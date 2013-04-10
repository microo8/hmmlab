#!/usr/bin/env python3
'''
This file is part of HMMLab.

HMMLab is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

HMMLab is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with HMMLab.  If not, see <http://www.gnu.org/licenses/>.
'''

import sys, os
import threading
from os.path import expanduser, join, exists
import configparser
from gi.repository import Gtk, Gdk
import cairo
try:
    from hmmlablib import libhmm
    from visual_win import VisualWindow
    import gtklib
    from model_win import ModelWindow
except ImportError:
    from hmmlab.hmmlablib import libhmm
    from hmmlab.visual_win import VisualWindow
    from hmmlab import gtklib
    from hmmlab.model_win import ModelWindow

class Select2DWindow(gtklib.ObjGetter):
    def __init__(self, modelset, main_window):
        path = join(os.path.dirname(os.path.abspath(__file__)), 'glade')
        gtklib.ObjGetter.__init__(self, join(path, 'select_2d.glade'), self.get_signals())
        self.modelset = modelset
        self.adjustment1.set_upper(self.modelset.streams_size)
        self.adjustment2.set_upper(self.modelset.streams_distribution[0])
        self.window.set_transient_for(main_window)
        self.window.show()

    def get_signals(self):
        signals = {"str_changed" : self.str_changed,
                "show_graph" : self.show_graph }
        return signals

    def destroy(self):
        self.window.destroy()

    def __del__(self):
        self.destroy()

    def str_changed(self, spin_button):
        self.adjustment2.set_upper(self.modelset.streams_distribution[int(self.adjustment1.get_value()-1)])

    def show_graph(self, button):
        stream_index = int(self.adjustment1.get_value() - 1)
        dim = int(self.adjustment2.get_value() - 1)
        self.modelset.gnuplot_2D(stream_index, dim)
        self.destroy()

class CanvasModel:
    def __init__(self, model, x, y, reset):
        self.model = model
        self.x = x
        self.y = y
        self.checked = False
        self.reset = reset

    def __repr__(self):
        return "%s, %s, %.2f, %.2f" % (self.model, self.model.name, self.x, self.y)

    def gaussians(self):
        if self.checked:
            self.model.select_gaussians()
        else:
            self.model.unselect_gaussians()
        self.reset()


class FilesTab(gtklib.ObjGetter):
    def __init__(self, modelset, main_window):
        path = join(os.path.dirname(os.path.abspath(__file__)), 'glade')
        gtklib.ObjGetter.__init__(self, join(path, 'files_tab.glade'), self.get_signals())
        self.modelset = modelset
        self.main_window = main_window

    def get_signals(self):
        signals = {"toggle_file" : self.toggle_file,
                   "destroy" : self.destroy,
                   "play" : self.play}
        return signals
    
    def destroy(self, widget, event):
        self.window.hide()
        return 1
    
    def toggle_file(self, widget, path):
        self.liststore[path][1] = not self.liststore[path][1]
        if self.liststore[path][1]:
            self.modelset.select_data(self.liststore[path][0])
        else:
            self.modelset.unselect_data(self.liststore[path][0])
        self.main_window.visual_win.refresh()
        self.load_data()

    def load_data(self):
        self.liststore.clear()
        self.liststore1.clear()
        for filename in self.modelset.files_data:
            fd = self.modelset.files_data[filename]
            self.liststore.append([filename, fd.selected, fd.word, fd.model.name if fd.model is not None else '', "%.6e" % (fd.maxprob if fd.model is not None else 0.0)])
        for model in self.modelset.drawarea_models:
            this_model_files = [self.modelset.files_data[x] for x in self.modelset.files_data if self.modelset.files_data[x].word != '' and self.modelset.files_data[x].word == model.name and self.modelset.files_data[x].model is not None]
            success_files = [x for x in this_model_files if x.word == x.model.name]
            if len(this_model_files) == 0:
                self.liststore1.append([model.name, "100%"])
            else:
                success = float(len(success_files)) / float(len(this_model_files)) * 100.0
                self.liststore1.append([model.name, "%.2f%%" % success])

    def play(self, treeview, path, col):
        threading.Thread(target=os.system, args=('aplay ' + self.liststore[path][0],)).start()
        
class MainWindow(gtklib.ObjGetter):
    '''Trieda hlavneho okna'''
    def __init__(self, modelset=None, modelset_path=None):
        '''Vytvori hlavne okno, nacita konfiguracny subor a nastavy velkost'''
        path = join(os.path.dirname(os.path.abspath(__file__)), 'glade')
        gtklib.ObjGetter.__init__(self, join(path, 'main.glade'), self.get_signals())
        self.config = configparser.ConfigParser()
        self.config.read(expanduser('~/.config/hmmlab.conf'))
        self.window.set_default_size(int(self.config['mainwindow']['width']), int(self.config['mainwindow']['height']))
        self.file_context_id = self.statusbar.get_context_id("Načítavanie súboru")
        self.MODEL_WIDTH = int(self.config['model']['width'])
        self.MODEL_HEIGHT = int(self.config['model']['height'])
        self.MODEL_CHECK = cairo.ImageSurface.create_from_png(join(path, 'check.png'))
        self.MODEL_UNCHECK = cairo.ImageSurface.create_from_png(join(path, 'uncheck.png'))
        self.mice = {'mouse_down':False, 'x': 0, 'y': 0, 'drag' : False}
        self.selection_rectangle = {'x1' : 0, 'y1' : 0, 'x2' : 0, 'y2' : 0} 
        
        self.modelset = modelset
        self.modelset_path = modelset_path
        self.data_table = FilesTab(modelset, self)
        self.modelset_modified = False
        self.visual_win = VisualWindow(self, self.modelset)
        if self.modelset is not None:
            self.fill_models()
            self.models_sidebar.set_sensitive(True)
            self.opened()
        self.edited = False

        self.models_canvas = []
        self.models_selected = []
        self.models_windows = []
        self.treeview1.drag_source_set(Gdk.ModifierType.BUTTON1_MASK, [], Gdk.DragAction.COPY)
        self.drawarea.drag_dest_set(Gtk.DestDefaults.ALL, [], Gdk.DragAction.COPY)
        self.treeview1.drag_source_add_text_targets()
        self.drawarea.drag_dest_add_text_targets()
        
        self.window.show()

    def get_signals(self):
        signals = {"about" : self.about,
                   "destroy" : self.destroy,
                   "drawarea_button_press_event_cb" : self.drawarea_button_press_event, 
                   "drawarea_button_release_event_cb" : self.drawarea_button_release_event,
                   "drawarea_draw_cb" : self.draw,
                   "drawarea_drag_data_received_cb" : self.drawarea_drag_data_received, 
                   "drawarea_motion_notify_event_cb" : self.drawarea_motion_notify_event,
                   "export_activate" : self.export_activate,
                   "gnuplot" : self.gnuplot,
                   "gnuplot3d" : self.gnuplot3D,
                   "load_data" : self.load_data,
                   "models_drag_get" : self.models_drag_get,
                   "open_activate" : self.open_activate,
                   "save_activate" : self.save_activate,
                   "save_as_activate" : self.save_as_activate,
                   "open_data_table" : self.open_data_table}
        return signals

    def destroy(self, widget = None, event=None):
        if widget is not None and self.modelset is not None and self.edited:
            dialog = Gtk.MessageDialog(self.window,
                                       0,
                                       Gtk.MessageType.QUESTION,
                                       Gtk.ButtonsType.YES_NO,
                                       "Uložiť zmeny do súboru '%s' pred zatvorením?" % self.modelset.name)
            dialog.format_secondary_text("Ak neuložíte, zmeny budú stratené.")
            response = dialog.run()
            if response == Gtk.ResponseType.YES:
                #bude treba otvorit save okno
                pass
            dialog.destroy()
        self.config['mainwindow']['width'] = str(self.window.get_allocated_width())
        self.config['mainwindow']['height'] = str(self.window.get_allocated_height())
        with open(expanduser('~/.config/hmmlab.conf'), 'w') as configfile:
            self.config.write(configfile)
        Gtk.main_quit()

    def draw(self, drawarea, cr):
        if self.modelset is not None:
            for model in self.models_canvas:
                if model in self.models_selected:
                    cr.set_source_rgba(0,80,0,0.5)
                else:
                    cr.set_source_rgba(0,0,0,0.6)
                #cr.rectangle(model.x, model.y, self.MODEL_WIDTH, self.MODEL_HEIGHT)
                #cr.fill()
                gtklib.cairo_rounded_rectangle(cr, model.x, model.y, self.MODEL_WIDTH, self.MODEL_HEIGHT, 1, self.MODEL_HEIGHT/10)
                cr.fill_preserve()
                cr.fill()
                cr.select_font_face("Sans", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
                cr.set_font_size(12.0)
                extents = cr.text_extents(model.model.name)
                x = model.x + self.MODEL_WIDTH / 2 - extents[2] / 2 - extents[0]
                y = model.y + self.MODEL_HEIGHT / 2 - extents[3] / 2 - extents[1]
                cr.move_to(x, y)
                cr.set_source_rgb(255,255,255)
                cr.show_text(model.model.name)
                if model.checked:
                    cr.set_source_surface(self.MODEL_CHECK, model.x+1, model.y+1)
                else:
                    cr.set_source_surface(self.MODEL_UNCHECK, model.x+1, model.y+1)
                cr.paint()
                
            if self.mice['mouse_down'] and not self.mice['drag']:
                cr.set_source_rgba(0,0,180,0.1)
                cr.rectangle(self.selection_rectangle['x1'],
                             self.selection_rectangle['y1'],
                             self.selection_rectangle['x2'] - self.selection_rectangle['x1'],
                             self.selection_rectangle['y2'] - self.selection_rectangle['y1'])
                cr.fill()
                cr.set_source_rgba(0,0,180,0.5)
                cr.rectangle(self.selection_rectangle['x1'],
                             self.selection_rectangle['y1'],
                             self.selection_rectangle['x2'] - self.selection_rectangle['x1'],
                             self.selection_rectangle['y2'] - self.selection_rectangle['y1'])
                cr.stroke()

    def opened(self):
        self.imagemenuitem3.set_sensitive(True)
        self.imagemenuitem4.set_sensitive(True)
        self.imagemenuitem11.set_sensitive(True)
        self.action1.set_sensitive(True)
        self.action2.set_sensitive(True)
        self.action3.set_sensitive(True)
        self.modelset_modified = False

    def save_activate(self, item):
        if self.modelset_modified:
            if self.modelset is None:
                self.save('xml')
            else:
                self.modelset.save(self.modelset_path, 'xml')

    def save_as_activate(self, item=None):
        self.save('xml')
        self.modelset_modified = False

    def export_activate(self, item=None):
        self.save('htk')

    def save(self, file_format):
        dialog = Gtk.FileChooserDialog("Vyberte súbor",
                                       self.window,
                                       Gtk.FileChooserAction.SAVE,
                                       (Gtk.STOCK_CANCEL,
                                       Gtk.ResponseType.CANCEL,
                                       Gtk.STOCK_OPEN,
                                       Gtk.ResponseType.OK))
        dialog.set_do_overwrite_confirmation(True)
        self.add_filters(dialog, file_format)
        response = dialog.run()
        if response == Gtk.ResponseType.OK:
            filename = dialog.get_filename()
            dialog.destroy()
            self.modelset.save(filename, file_format)
        else:
            dialog.destroy()

    def load_data(self, item):
        dialog = Gtk.FileChooserDialog("Vyberte súbor",
                                       self.window,
                                       Gtk.FileChooserAction.OPEN,
                                       (Gtk.STOCK_CANCEL,
                                       Gtk.ResponseType.CANCEL,
                                       Gtk.STOCK_OPEN,
                                       Gtk.ResponseType.OK))
        filter_wav = Gtk.FileFilter()
        filter_wav.set_name("WAV súbor")
        filter_wav.add_pattern('*.wav')
        dialog.add_filter(filter_wav)
        dialog.set_select_multiple(True)
        response = dialog.run()
        if response == Gtk.ResponseType.OK:
            filenames = dialog.get_filenames()
            self.modelset.load_data(filenames)
            self.visual_win.refresh()
        dialog.destroy()

    def open_data_table(self, item):
        self.data_table.modelset = self.modelset
        self.data_table.load_data()
        self.data_table.window.show()

    def open_activate(self, item):
        dialog = Gtk.FileChooserDialog("Vyberte súbor",
                                        self.window,
                                        Gtk.FileChooserAction.OPEN,
                                        (Gtk.STOCK_CANCEL,
                                        Gtk.ResponseType.CANCEL,
                                        Gtk.STOCK_OPEN,
                                        Gtk.ResponseType.OK))
        self.add_filters(dialog)
        response = dialog.run()
        if response == Gtk.ResponseType.OK:
            filename = dialog.get_filename()
            file_type = filename[filename.index('.')+1:]
            dialog.destroy()
#vyriesit aby sa to ukazalo v statusbare
            self.statusbar.push(self.file_context_id,'Načítavam súbor '+filename)
            if self.modelset is not None:
                self.modelset.destroy()
            self.modelset = libhmm.ModelSet(filename, file_type)
            self.visual_win.set_modelset(self.modelset)
            self.statusbar.pop(self.file_context_id)
            self.window.set_title('HMMLab - '+self.modelset.name)
            self.drawarea.queue_draw()
            self.fill_models()
            self.models_sidebar.set_sensitive(True)
            self.opened()
        else:
            dialog.destroy()

    def add_filters(self, dialog, file_format=None):
        if file_format == 'htk' or file_format is None:
            filter_hmm = Gtk.FileFilter()
            filter_hmm.set_name("HTK HMM súbor")
            filter_hmm.add_pattern('*.hmm')
            dialog.add_filter(filter_hmm)

        if file_format == 'xml' or file_format is None:
            filter_xml = Gtk.FileFilter()
            filter_xml.set_name("XML súbory")
            filter_xml.add_mime_type("text/xml")
            dialog.add_filter(filter_xml)

        filter_any = Gtk.FileFilter()
        filter_any.set_name("Všetky súbory")
        filter_any.add_pattern("*")
        dialog.add_filter(filter_any)

    def fill_models(self):
        if self.modelset is not None:
            self.models_store.clear()
            for model in self.modelset.models:
                self.models_store.append([model.name])

    def models_drag_get(self, widget, drag_context, data, info, time):
        if self.modelset is not None:
            model, it = self.treeview1.get_selection().get_selected()
            item = self.models_store.get_value(it, 0)
            data.set_text(item, -1)

    def drawarea_drag_data_received(self, widget, drag_context, x, y, data, info, time):
        if self.modelset is not None:
            model_name = data.get_text()
            model = self.modelset.get_model(model_name)
            self.models_canvas.append(CanvasModel(model, x, y, self.modelset.reset_pos_gauss))
            self.modelset.drawarea_models_append(model)
            self.drawarea.queue_draw()
            self.modelset_modified = True

    def drawarea_button_press_event(self, eb, event):
        if self.modelset is not None:
            model = self.get_model_at(event.x, event.y)
            if model is None:
                self.models_selected = []
                self.mice['drag'] = False
                self.selection_rectangle['x1'] = event.x
                self.selection_rectangle['y1'] = event.y
                self.selection_rectangle['x2'] = event.x
                self.selection_rectangle['y2'] = event.y
            else:
                if event.type == Gdk.EventType._2BUTTON_PRESS:
                    if model.model.name not in [mw.model.name for mw in self.models_windows]:
                        mw = ModelWindow(model, self)
                        self.models_windows.append(mw)
                    else:
                        for mw in self.models_windows:
                            if model.model.name == mw.model.name:
                                mw.window.present()
                                break
                if model not in self.models_selected:
                    self.models_selected = [model]
                if model.x <= event.x <= model.x + self.MODEL_CHECK.get_width():
                    if model.y <= event.y <= model.y + self.MODEL_CHECK.get_height():
                        model.checked = not model.checked
                        model.gaussians()                        
                        for mw in self.models_windows:
                            if model.model.name == mw.model.name:
                                mw.fill_states_table()
                        self.visual_win.refresh()
                self.mice['drag'] = True
            self.mice['mouse_down'] = True
            self.mice['x'] = event.x
            self.mice['y'] = event.y
            self.drawarea.queue_draw()

    def drawarea_button_release_event(self, eb, event):
        if self.modelset is not None:
            self.mice['mouse_down'] = False
            self.mice['x'] = event.x
            self.mice['y'] = event.y
            if not self.mice['drag']:
                self.update_selected_models()
            self.mice['drag'] = False
            for k in self.selection_rectangle:
                self.selection_rectangle[k] = -1
            self.drawarea.queue_draw()

    def drawarea_motion_notify_event(self, eb, event):
        if self.modelset is not None:
            if self.mice['mouse_down']:
                alloc = self.drawarea.get_allocation()
                if self.mice['drag']:
                    for model in self.models_selected:
                        model.x = max(0, min(alloc.width - self.MODEL_WIDTH, model.x - self.mice['x'] + event.x))
                        model.y = max(0, min(alloc.height - self.MODEL_HEIGHT, model.y - self.mice['y'] + event.y))
                    if len(self.models_selected) > 0:
                        self.drawarea.queue_draw()
                elif len(self.models_selected) == 0:
                        self.selection_rectangle['x2'] = max(0, min(alloc.width - self.MODEL_WIDTH,event.x))
                        self.selection_rectangle['y2'] = max(0, min(alloc.height - self.MODEL_HEIGHT, event.y))
                        self.drawarea.queue_draw()
            self.mice['x'] = event.x
            self.mice['y'] = event.y
        
    def update_selected_models(self):
        self.models_selected = []
        for model in self.models_canvas:
            if self.selection_rectangle['x1'] > self.selection_rectangle['x2']:
                if self.selection_rectangle['x1'] >= (model.x+self.MODEL_WIDTH) and self.selection_rectangle['x2'] <= model.x:
                    self.models_selected.append(model)
            else:
                if self.selection_rectangle['x2'] >= (model.x+self.MODEL_WIDTH) and self.selection_rectangle['x1'] <= model.x:
                    self.models_selected.append(model)

    def get_model_at(self, x, y):
        for model in self.models_canvas:
            if model.x <= x <= model.x + self.MODEL_WIDTH and model.y <= y <= model.y + self.MODEL_HEIGHT:
                return model
        return None

    def about(self, item):
        resp = self.aboutdialog.run()
        if resp == Gtk.ResponseType.DELETE_EVENT or resp == Gtk.ResponseType.CANCEL:
            self.aboutdialog.hide()

    def gnuplot(self, item):
        Select2DWindow(self.modelset, self.window)

    def gnuplot3D(self, item):
        SelectDraw2DWindow(self.modelset, self)
    
    def open_window_from_gauss(self, gauss):
        for model in self.modelset.get_models_with_gaussian(gauss):
            mw = None
            if model.name not in [mw.model.name for mw in self.models_windows]:
                mw = ModelWindow(model, self)
                self.models_windows.append(mw)
            else:
                for mwin in self.models_windows:
                    if model.name == mwin.model.name:
                        mwin.window.present()
                        mw = mwin
                        break

            if mw is not None:
                for i, state in enumerate(model.states):
                    if state.has_gaussian(gauss):
                        selection = mw.states_view.get_selection()
                        if selection is not None:
                            selection.select_path(i)
                            mw.load_state(state)
                            for j, data in enumerate(mw.gaussians_store):
                                if data[2] == gauss.name:
                                    selection = mw.gaussians_view.get_selection()
                                    if selection is not None:
                                        selection.select_path(j)


def run():
    if len(sys.argv) > 1:
        if exists(sys.argv[1]):
            file_type = sys.argv[1][sys.argv[1].index('.')+1:]
            modelset = libhmm.ModelSet(sys.argv[1], file_type)
            if len(sys.argv) > 3 and sys.argv[2] == '-w':
                filenames = [filename for filename in sys.argv[3:] if exists(filename)]
                modelset.load_data(filenames)
            MainWindow(modelset)
            Gtk.main()
        else:
            print("Súbor", sys.argv[1], "neexistuje")
    else:
        MainWindow()
        Gtk.main()

if __name__ == '__main__':
    run()
