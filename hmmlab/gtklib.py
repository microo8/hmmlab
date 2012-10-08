#    This file is part of HMMLab.
#
#    HMMLab is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    HMMLab is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with HMMLab.  If not, see <http://www.gnu.org/licenses/>.

from gi.repository import Gtk

class ObjGetter:
    def __init__(self, path, signals):
        self.builder = Gtk.Builder()
        self.builder.add_from_file(path)
        self.builder.connect_signals(signals)

    def __getattr__(self, name):
        return self.builder.get_object(name)
