#!/usr/bin/env python3

<<<<<<< HEAD
import sys, os
=======
<<<<<<< HEAD
import sys
=======
import sys, pickle, os
>>>>>>> new_branch_name
>>>>>>> parent of af2eb34... no uz prestan blbnut
from os.path import expanduser, join, exists
import configparser
from gi.repository import Gtk, Gdk
import hmmlablib
import gtklib
import sqlalchemy as sql

class CanvasModel:
    def __init__(self, model, x, y):
        self.model = model
        self.x = x
        self.y = y
        self.z = 0

    def __repr__(self):
        return "%s, %.2f, %.2f" % (self.model.name, self.x, self.y)

class MainWindow(gtklib.ObjGetter):
    '''Trieda hlavneho okna'''
    def __init__(self, modelset=None):
        '''Vytvori hlavne okno, nacita konfiguracny subor a nastavy velkost'''
        gtklib.ObjGetter.__init__(self, join('glade', 'main.glade'), self.get_signals())
        self.config = configparser.ConfigParser()
        self.config.read(expanduser('~/.config/hmmlab.conf'))
        self.window.set_default_size(int(self.config['mainwindow']['width']), int(self.config['mainwindow']['height']))
        self.file_context_id = self.statusbar.get_context_id("Načítavanie súboru")
        self.MODEL_WIDTH = int(self.config['model']['width'])
        self.MODEL_HEIGHT = int(self.config['model']['height'])
        self.mice = {'mouse_down':False, 'x': 0, 'y': 0, 'drag' : False}
        self.selection_rectangle = {'x1' : 0, 'y1' : 0, 'x2' : 0, 'y2' : 0} 
        
        self.modelset = modelset
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
        signals = {"destroy" : self.destroy,
                   "drawarea_draw_cb" : self.draw,
                   "open_activate" : self.open_activate,
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
                                       Gtk.ButtonsType.YES_NO, "Uložiť zmeny do súboru '%s' pred zatvorením?" % self.modelset.name)
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
                cr.rectangle(model.x, model.y, self.MODEL_WIDTH, self.MODEL_HEIGHT)
                if model in self.models_selected:
                    cr.set_source_rgba(0,80,0,0.5)
                else:
                    cr.set_source_rgba(0,0,0,0.5)
                cr.fill()
            if self.mice['mouse_down'] and not self.mice['drag']:
                cr.rectangle(self.selection_rectangle['x1'],
                             self.selection_rectangle['y1'],
                             self.selection_rectangle['x2'] - self.selection_rectangle['x1'],
                             self.selection_rectangle['y2'] - self.selection_rectangle['y1'])
                cr.set_source_rgba(0,0,180,0.1)
                cr.fill()
                cr.rectangle(self.selection_rectangle['x1'],
                             self.selection_rectangle['y1'],
                             self.selection_rectangle['x2'] - self.selection_rectangle['x1'],
                             self.selection_rectangle['y2'] - self.selection_rectangle['y1'])
                cr.set_source_rgba(0,0,180,0.5)
                cr.stroke()

    def open_activate(self, item):
        dialog = Gtk.FileChooserDialog("Vyberte súbor", self.window,
            Gtk.FileChooserAction.OPEN,
            (Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL,
             Gtk.STOCK_OPEN, Gtk.ResponseType.OK))

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
            self.modelset = hmmlablib.ModelSet(filename, file_type)
            self.statusbar.pop(self.file_context_id)
            self.window.set_title('HMMLab - '+self.modelset.name)
            self.drawarea.queue_draw()
            self.fill_models()
            self.models_sidebar.set_sensitive(True)
        else:
            dialog.destroy()

    def add_filters(self, dialog):
        filter_hmm = Gtk.FileFilter()
        filter_hmm.set_name("HTK HMM súbor")
        filter_hmm.add_pattern('*.hmm')
        dialog.add_filter(filter_hmm)

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
                if self.mice['drag']:
                    for model in self.models_selected:
                        model.x -= self.mice['x'] - event.x
                        model.y -= self.mice['y'] - event.y
                    if len(self.models_selected) > 0:
                        self.drawarea.queue_draw()
                elif len(self.models_selected) == 0:
                        self.selection_rectangle['x2'] = event.x
                        self.selection_rectangle['y2'] = event.y
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
        

if __name__ == '__main__':
    if len(sys.argv) > 1:
        if exists(sys.argv[1]):
            file_type = sys.argv[1][sys.argv[1].index('.')+1:]
            modelset = hmmlablib.ModelSet(sys.argv[1], file_type)
            MainWindow(modelset)
            Gtk.main()
        else:
            print("Súbor", sys.argv[1], "neexistuje")
    else:
        MainWindow()
        Gtk.main()
