# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 08:06:37 2020

@author: Connor
"""

from tkinter import *
from tkinter.filedialog import askopenfilenames, asksaveasfilename
import hyperspy.api as hs
import numpy as np
from Load_File import load_file

def tk_param():
    root = Tk()
    root.title("Holography Reconstruction Parameters")
    
    def set1(sb_m):
        #global man1
        #global frame3
        #global sb_position
        if sb_m==False:
            man1 = Entry(frame3, textvariable=sb_position, state=DISABLED, width=12)
        else:
            man1 = Entry(frame3, textvariable=sb_position, state=NORMAL, width=12)
        man1.grid(row=2, column=1, sticky=W)
        
    def set2(ap_m):
        #global man2
        #global frame4
        #global ap_size
        if ap_m==False:
            man2 = Entry(frame4, textvariable=ap_size, state=DISABLED, width=5)
            ph_weak = Radiobutton(frame4, text='Weak Phase Object', state=NORMAL, 
                                  padx=10, variable=ph_s, value=False)
            ph_strong = Radiobutton(frame4, text='Strong Phase Object', 
                                    state= NORMAL, padx=10, variable=ph_s, 
                                    value=True)
        else:
            man2 = Entry(frame4, textvariable=ap_size, state=NORMAL, width=5)
            ph_weak = Radiobutton(frame4, text='Weak Phase Object', state=DISABLED, 
                                  padx=10, variable=ph_s, value=False)
            ph_strong = Radiobutton(frame4, text='Strong Phase Object', 
                                    state= DISABLED, padx=10, variable=ph_s, 
                                    value=True)
        man2.grid(row=4, column=1, sticky=W)
        ph_strong.grid(row=3, column=0, columnspan=2, sticky=W)
        ph_weak.grid(row=2, column=0, columnspan=2, sticky=W)
    
    def set3(butterworth):
        #global man3
        #global frame5
        #global butter_var
        if butterworth==False:
            man3 = Entry(frame5, textvariable=butter_var, state=DISABLED, width=5)
        else:
            man3 = Entry(frame5, textvariable=butter_var, state=NORMAL, width=5)
        man3.grid(row=3, column=1)
        
    def set4(manual):
        #global man4
        #global frame6
        #global rec_var
        if manual==False:
            man4 = Entry(frame6, textvariable=rec_var, state=DISABLED, width=5)
        else:
            man4 = Entry(frame6, textvariable=rec_var, state=NORMAL, width=5)
        man4.grid(row=3, column=1, sticky=W)
    
    def ref_ask():
        nonlocal ref_files
        #global ref_label
        ref_files = askopenfilenames(title="Choose reference files", filetypes=[("", "*.dm3")])
        #ref_label.grid_forget()
        Label(frame1, text="reference files selected", fg='green').grid(row=3, column=0, columnspan=2, sticky=W)
        
    def obj_ask():
        nonlocal obj_files
        #global obj_label
        obj_files = askopenfilenames(title="Choose object files", filetypes=[("", "*.dm3")])
        #obj_label.grid_forget()
        Label(frame1, text="object files selected", fg='green').grid(row=1, column=0, columnspan=2, sticky=W)
        
    def out_ask():
        nonlocal out_file
        out_file = asksaveasfilename(title="Save output as:", defaultextension=".hdf5", filetypes=[("", "*.hdf5")])
        Label(frame7, text=out_file, fg='green').grid(row=1, column=0, sticky=W, columnspan=3)
        
    def quit_tk():
        root.destroy()
        sys.exit("Reconstruction canceled.")
        
    def run_tk():
        root.destroy()
        
    def get_mtf():
        return 1
    
    def get_sample(file):
        #convert from pixel radius to frequency cutoff
        sample = load_file(file).attrs["dim_scale"]
        return sample[0]
    
    # Initialize variables
    align_var = BooleanVar()
    fringe_var = BooleanVar()
    defocus_var = BooleanVar()
    res_var = BooleanVar()
    fresnel_var = BooleanVar()
    mtf_var = BooleanVar()
    sb_m = BooleanVar()
    sb_position = StringVar()
    ap_m = BooleanVar()
    ap_size = DoubleVar()
    ph_s = BooleanVar()
    butter_var = DoubleVar()
    ap_type = StringVar()
    re_size = DoubleVar()
    rec_var = IntVar()
    man4_var = IntVar()
    ref_files = ()
    obj_files = ()
    out_file = ""
    
    # Set reference and object directories
    frame1 = Frame(root, padx=10, pady=5)
    frame1.grid(row=0, column=0, columnspan=2, sticky=W)
    
    ref_button = Button(frame1, text='Choose reference files',command=ref_ask)
    obj_button = Button(frame1, text='Choose object files',command=obj_ask)
   
    obj_button.grid(row=0, column=0, sticky=W)
    Label(frame1, text=" ").grid(row=1, column=0, sticky=W)
    ref_button.grid(row=2, column=0, sticky=W)
    Label(frame1, text=" ").grid(row=3, column=0, sticky=W)
    
    # Determine which corrections to apply
    frame2 = Frame(root, padx=10, pady=5)#, bd=1, relief=SUNKEN)
    frame2.grid(row=1, column=0, columnspan=2, sticky=W)
    
    Label(frame2, text="Corrections to Apply").grid(row=0, column=0, columnspan=2, 
                                                    sticky=W)
    align_check = Checkbutton(frame2, text='Sample drift', variable=align_var)
    fringe_check = Checkbutton(frame2, text='Fringe drift', variable=fringe_var)
    defocus_check = Checkbutton(frame2, text='Defocus', variable=defocus_var)
    res_check = Checkbutton(frame2, text='High resolution', variable=res_var)
    fresnel_check = Checkbutton(frame2, text='Biprism Fresnel fringes', 
                                variable=fresnel_var)
    mtf_check = Checkbutton(frame2, text='MTF', variable=mtf_var)
    
    align_check.grid(row=1, column=0, sticky=W)
    fringe_check.grid(row=2, column=0, sticky=W)
    defocus_check.grid(row=3, column=0, sticky=W)
    res_check.grid(row=1, column=1, sticky=W)
    fresnel_check.grid(row=2, column=1, sticky=W)
    mtf_check.grid(row=3, column=1, sticky=W)
    
    
    # Sideband position
    frame3 = Frame(root, padx=10, pady=5)
    frame3.grid(row=2, column=0, sticky=N)
    
    Label(frame3, text="Sideband Position (FFT px)").grid(row=0, column=0, 
                                                              columnspan=2)
    sb_man = Radiobutton(frame3, text='Manual [x, y]', variable=sb_m, value=True, 
                         command=lambda: set1(True))
    sb_auto = Radiobutton(frame3, text='Auto', variable=sb_m, value=False, 
                          command=lambda: set1(False))
    man1 = Entry(frame3, textvariable=sb_position, state=DISABLED, width=12)
    
    sb_auto.grid(row=1, column=0, sticky=W)
    sb_man.grid(row=2, column=0, sticky=W)
    man1.grid(row=2, column=1, sticky=W)
    sb_m.set(False)
    
    # Aperture Cutoff
    frame4 = Frame(root, padx=10, pady=5)
    frame4.grid(row=2, column=1)
    
    Label(frame4, text="Aperture Cutoff (1/nm)").grid(row=0, column=0, columnspan=2, 
                                                     sticky=W)
    ap_man = Radiobutton(frame4, text='Manual', variable=ap_m, value=True, 
                         command=lambda: set2(True))
    ap_auto = Radiobutton(frame4, text='Auto', variable=ap_m, value=False, 
                          command=lambda: set2(False))
    ph_weak = Radiobutton(frame4, text='Weak Phase Object', padx=10, variable=ph_s, 
                             value=False)
    ph_strong = Radiobutton(frame4, text='Strong Phase Object', padx=10, variable=ph_s, 
                               value=True)
    man2 = Entry(frame4, textvariable=ap_size, state=DISABLED, width=5)
    
    ph_strong.grid(row=3, column=0, columnspan=2, sticky=W)
    ph_weak.grid(row=2, column=0, columnspan=2, sticky=W)
    ap_man.grid(row=4, column=0, sticky=W)
    ap_auto.grid(row=1, column=0, sticky=W)
    man2.grid(row=4, column=1, sticky=W)
    ap_m.set(False)
    ph_s.set(False)
    
    # Aperture Type
    frame5 = Frame(root, padx=10, pady=5)
    frame5.grid(row=3, column=0, sticky=N+W)
    
    Label(frame5, text="Aperture Type").grid(row=0, column=0, sticky=W)
    ap_E = Radiobutton(frame5, text='Hard Edge', variable=ap_type, value="EDGE", 
                         command=lambda: set3(False))
    ap_G = Radiobutton(frame5, text='Gaussian', variable=ap_type, value="GAUSSIAN", 
                          command=lambda: set3(False))
    ap_B = Radiobutton(frame5, text='Butterworth', variable=ap_type, 
                       value="BUTTERWORTH", command=lambda: set3(True))
    man3 = Entry(frame5, textvariable=butter_var, state=DISABLED, width=5)
    ap_type.set("EDGE")
    
    ap_E.grid(row=1, column=0, sticky=W)
    ap_G.grid(row=2, column=0, sticky=W)
    ap_B.grid(row=3, column=0, sticky=W)
    man3.grid(row=3, column=1, sticky=W)
    
    # Reconstruction Size
    frame6 = Frame(root, padx=10, pady=5)
    frame6.grid(row=3, column=1, sticky=W)
    
    Label(frame6, text="Reconstructed Size (px)").grid(row=0, column=0, 
         columnspan=2, sticky=W)
    re_2 = Radiobutton(frame6, text='1/2 resize', variable=re_size, value=0.5, 
                         command=lambda: set4(False))
    re_4 = Radiobutton(frame6, text='1/4 resize', variable=re_size, value=0.25, 
                          command=lambda: set4(False))
    re_man = Radiobutton(frame6, text='Manual', variable=re_size, 
                       value=man4_var, command=lambda: set4(True))
    man4 = Entry(frame6, textvariable=man4_var, state=DISABLED, width=5)
    re_size.set(0.5)
    
    re_2.grid(row=1, column=0, columnspan=2, sticky=W)
    re_4.grid(row=2, column=0, columnspan=2, sticky=W)
    re_man.grid(row=3, column=0, sticky=W)
    man4.grid(row=3, column=1, sticky=W)
    
    # Save file, Quit, Run buttons
    frame7 = Frame(root, padx=10, pady=5)
    frame7.grid(row=4, column=0, columnspan=2, sticky=E+W)
    
    out_button = Button(frame7, text="Save output as:", command=out_ask)
    out_button.grid(row=0, column=0, sticky=W)
    Label(frame7, text="").grid(row=1, column=0, sticky=W)
    quit_button = Button(frame7, text="QUIT", command=quit_tk, bg='red', fg='white', padx=5)
    run_button = Button(frame7, text="START", command=run_tk, bg='green', fg='white', padx=5)
    frame7.grid_columnconfigure(1, weight=1)
    quit_button.grid(row=0, column=1, sticky=E, padx=5)
    run_button.grid(row=0, column=2, sticky=E, padx=5)
    
    root.mainloop()
    
    """Store parameters from menu in dictionary
    Use HyperSpy for some auto settings
    Not yet implemented: MTF input"""
    
    ref_holo = hs.load(ref_files[0], signal_type='hologram')
    holo_size = (ref_holo.data.shape[0], ref_holo.data.shape[1])
    
    if sb_m.get()==False:
        #get sideband position and size in HyperSpy format, convert later
        sb_position = ref_holo.estimate_sideband_position()
        sb_position = sb_position.data
        
        for i in range(2):
            if sb_position[i] > holo_size[i]/2:
                sb_position[i] -= holo_size[i]/2
            else:
                sb_position[i] += holo_size[i]/2
                
        sb_position = np.flip(sb_position)
        sb_position = sb_position.astype(int)
        sb_position = [sb_position[0].item(), sb_position[1].item()]
    else:
        sb_position = sb_position.get().replace('[','').replace(']','').split(',')
        sb_position = [int(sb_position[0]), int(sb_position[1])]
        
    carry_freq = np.array([((sb_position[0] - holo_size[1]//2)/holo_size[1]), 
                  ((sb_position[1] - holo_size[0] // 2) / holo_size[0])])
    carry_freq = carry_freq[...] / get_sample(ref_files[0]) #get carry frequency
    carry_freq = np.sqrt(carry_freq[0] ** 2 + carry_freq[1] ** 2)
        
    if ap_m.get()==False:
        if ph_s.get()==True:
            ap_size = (1/3) * carry_freq
        else:
            ap_size = 0.5 * carry_freq   
    else:
        ap_size = ap_size.get()    
           
    if ap_type.get() == "BUTTERWORTH":
        ap_type = [ap_type.get(), butter_var.get()]
    else:
        ap_type = ap_type.get()
        
    if holo_size[1] < holo_size[0]:
        re_size = int(re_size.get() * holo_size[1])
    else:
        re_size = int(re_size.get() * holo_size[0])
    
    param = {'object_names': obj_files,
     'object_first': 0,
     'object_last': 0,
     'empty_names': ref_files,
     'empty_first': 0,
     'empty_last': 0,
     'defocus_first': 0,
     'empty_size': re_size,
     'object_size': re_size,
     'sideband_pos': sb_position,
     'output_name': out_file,
     'filter_func': ap_type,
     'cut_off': ap_size,
     'adjust_defocus': res_var.get(),
     'adjust_shift': res_var.get(),
     'adjust_tilt': fringe_var.get(),
     'bpf_correct': fresnel_var.get()}
    
    if mtf_var.get()==1:
        param.add('mtf', get_mtf())

    return param

if __name__=='__main__':
    print(tk_param())