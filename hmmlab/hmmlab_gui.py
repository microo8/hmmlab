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
from os.path import expanduser, join, exists
import configparser
from gi.repository import Gtk, Gdk
import cairo
try:
    from hmmlablib import libhmm
    from visual_win import VisualWindow
    import gtklib
except ImportError:
    from hmmlab.hmmlablib import libhmm
    from hmmlab.visual_win import VisualWindow
    from hmmlab import gtklib

class CanvasModel:
    def __init__(self, model, x, y):
        self.model = model
        self.x = x
        self.y = y
        self.checked = False

    def __repr__(self):
        return "%s, %.2f, %.2f" % (self.model.name, self.x, self.y)


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
        self.modelset_modified = False
        self.visual_win = VisualWindow(self.modelset)
        if self.modelset is not None:
            self.fill_models()
            self.models_sidebar.set_sensitive(True)
        self.edited = False

        self.models_canvas = []
        self.models_selected = []
        self.treeview1.drag_source_set(Gdk.ModifierType.BUTTON1_MASK, [], Gdk.DragAction.COPY)
        self.drawarea.drag_dest_set(Gtk.DestDefaults.ALL, [], Gdk.DragAction.COPY)
        self.treeview1.drag_source_add_text_targets()
        self.drawarea.drag_dest_add_text_targets()
        
        self.window.show()

    def get_signals(self):
        signals = {"about" : self.about,
                   "destroy" : self.destroy,
                   "drawarea_draw_cb" : self.draw,
                   "open_activate" : self.open_activate,
                   "save_activate" : self.save_activate,
                   "save_as_activate" : self.save_as_activate,
                   "export_activate" : self.export_activate,
                   "models_drag_get" : self.models_drag_get,
                   "drawarea_drag_data_received_cb" : self.drawarea_drag_data_received, 
                   "drawarea_button_press_event_cb" : self.drawarea_button_press_event, 
                   "drawarea_button_release_event_cb" : self.drawarea_button_release_event,
                   "drawarea_motion_notify_event_cb" : self.drawarea_motion_notify_event}
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
            self.models_canvas.append(CanvasModel(self.modelset.get_model(model_name), x, y))
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
                if model not in self.models_selected:
                    self.models_selected = [model]
                if model.x <= event.x <= model.x + self.MODEL_CHECK.get_width():
                    if model.y <= event.y <= model.y + self.MODEL_CHECK.get_height():
                        model.checked = not model.checked
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
        
def run():
    if len(sys.argv) > 1:
        if exists(sys.argv[1]):
            file_type = sys.argv[1][sys.argv[1].index('.')+1:]
            modelset = libhmm.ModelSet(sys.argv[1], file_type)
            MainWindow(modelset)
            Gtk.main()
        else:
            print("Súbor", sys.argv[1], "neexistuje")
    else:
        MainWindow()
        Gtk.main()

if __name__ == '__main__':
    run()
