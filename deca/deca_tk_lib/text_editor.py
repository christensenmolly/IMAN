#Taken from http://my-python3-code.blogspot.ru/2012/08/basic-tkinter-text-editor-online-example.html
#import tkinter.filedialog
import Tkinter
from Tkinter import *
from ScrolledText import *
import tkFileDialog
import tkMessageBox
    
class App:

    def new_command(self):
            # Clear the text
            self.textPad.delete(0.0, END)

    def open_command(self):
            file = tkFileDialog.askopenfile(parent=self.root, mode='rb', title='Select a file')
            if file != None:
	      contents = file.read()
	      self.textPad.insert('1.0',contents)
	      file.close()
    def save_command(self):
	    file = tkFileDialog.asksaveasfile(mode='w')
	    if file != None:
	      data = self.textPad.get('1.0', END+'-1c')
	      file.write(data)
	      file.close()
	      
    def exit_command():
	    if tkMessageBox.askokcancel("Quit", "Do you really want to quit?"):
	      self.root(destroy)

    def __init__(self,inp_file=None):
            # Set up the screen, the title, and the size.
            self.root = Tkinter.Tk(className="Text editor")
            self.textPad = ScrolledText(self.root, width=100, height=80)
      
	    menu = Menu(self.root)
	    self.root.config(menu=menu)
	    filemenu = Menu(menu)
	    menu.add_cascade(label="File", menu=filemenu)
	    filemenu.add_command(label="New", command=self.new_command)
	    filemenu.add_command(label="Open", command=self.open_command)
	    filemenu.add_command(label="Save", command=self.save_command)
	    filemenu.add_separator()
	    filemenu.add_command(label="Exit", command=self.exit_command)

	    
	    self.textPad.pack()
	    self.root.mainloop()

