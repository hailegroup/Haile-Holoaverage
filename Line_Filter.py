# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 13:42:58 2020

@author: Connor
"""
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.path import Path
import matplotlib.patches as patches
from matplotlib.widgets import Slider, Button

class LineFilter:
    
    def __init__(self, image, roi=None):
        self.fig, self.ax = plt.subplots()
        if roi is None:
            roi = [0, image.shape[1], 0, image.shape[0]]
        x1, x2, y1, y2 = roi[0], roi[1], roi[2], roi[3]
        self.image_full = image
        self.image_crop = image[x1:x2, y1:y2]
        self.roi = roi
        plt.xlim([x1, x2])
        plt.ylim([y2, y1])
        plt.subplots_adjust(bottom=0.1, left=0.1, top=0.9, right=0.9)
        self.ax.tick_params(axis='both', which='both', bottom=False, left=False,
                       labelbottom=False, labelleft=False)
        self.ax.imshow(self.image_crop, extent=[x1, x2, y2, y1])
        self.num = 0
        self.x = [0, 0]
        self.y = [0, 0]
        self.pts = np.ones(image.shape, dtype=bool)
        self.cid1 = self.fig.canvas.mpl_connect('button_press_event', self.onclick)
        self.cid2 = self.fig.canvas.mpl_connect('button_release_event', self.mouse_release)
        self.cid3 = self.fig.canvas.mpl_connect('motion_notify_event', self.mouse_move)
        self.cid4 = self.fig.canvas.mpl_connect('pick_event', self.onpick)
        self.cid5 = self.fig.canvas.mpl_connect('key_press_event', self.del_key)
        ax2 = plt.axes([0.2, 0.05, 0.4, 0.03], facecolor='lightgoldenrodyellow')
        self.lw = Slider(ax2, 'Line Width (px)', 1, 20, valinit=5, valfmt='%d', 
                         closedmin=True, closedmax=True, dragging=True, valstep=1, 
                         orientation='horizontal')
        axres = plt.axes([0.7, 0.04, 0.1, 0.05])
        self.bres = Button(axres, 'Reset', image=None)
        self.bres.on_clicked(self.reset)
        axacc = plt.axes([0.82, 0.04, 0.1, 0.05])
        self.bacc = Button(axacc, 'Accept', image=None)
        self.bacc.on_clicked(self.accept)
        plt.show()

    def accept(self, event):
        plt.close(self.fig)
    
    def reset(self, event):
        self.ax.clear()
        roi = self.roi
        x1, x2, y1, y2 = roi[0], roi[1], roi[2], roi[3]
        self.ax.imshow(self.image_crop, extent=[x1, x2, y2, y1])
        self.num = 0
        self.fig.canvas.draw_idle()
    
    def draw_line(self):
        self.line.set_data(self.x, self.y)
        self.fig.canvas.draw_idle()
    
    def mouse_move(self, event):
        if self.num == 1:
            self.x[1] = event.xdata
            self.y[1] = event.ydata
            self.draw_line()
    
    def mouse_release(self, event):
        if self.num == 1:
            self.num = 2                
    
    def onclick(self, event):
        if event.button == 1:
            # Check if mouse click is inside image
            if self.roi[0] < event.xdata < self.roi[1] and self.roi[2] < event.ydata < self.roi[3]:
                # Draw line
                if self.num==0:
                    self.num = 1
                    self.x = [event.xdata, event.xdata]
                    self.y = [event.ydata, event.ydata]
                    self.line, = self.ax.plot(self.x, self.y, 'r', picker=5)
                else:
                    pass
            else:
                pass
        else:
            pass # Do nothing
    
    def onpick(self, event):    
        """
        Handles the pick event - if an object has been picked, store a
        reference to it.  We do this by simply adding a reference to it
        named 'stored_pick' to the axes object.  Note that in python we
        can dynamically add an attribute variable (stored_pick) to an 
        existing object - even one that is produced by a library as in this
        case
        """
        this_artist = event.artist #the picked object is available as event.artist
        self.ax.picked_object = this_artist
    
    def del_key(self, event):
        """
        Function to be bound to the key press event
        If the key pressed is delete and there is a picked object,
        remove that object from the canvas
        """
        ax = self.ax    
        if event.key == u'delete':
            if ax.picked_object:
                self.num = 0
                ax.picked_object.remove()
                ax.picked_object = None
                ax.figure.canvas.draw()
            
        if event.key == u'enter':
            if self.num==2:
                ax.lines.pop(0)
                self.parallel_lines(self.lw.val)
                patch = patches.PathPatch(self.path, facecolor='k', lw=1)
                self.ax.add_patch(patch)
                self.fig.canvas.draw_idle()
                self.num = 0
        
        if event.key == 'q':
            self.close_connections()
                
    def close_connections(self):
        self.fig.canvas.mpl_disconnect(self.cid1)
        self.fig.canvas.mpl_disconnect(self.cid2)
        self.fig.canvas.mpl_disconnect(self.cid3)
        self.fig.canvas.mpl_disconnect(self.cid4)
        self.fig.canvas.mpl_disconnect(self.cid5)
                
    def parallel_lines(self, t):
        """input line coordinates and thickness
        return (rotated) rectangular path"""
        x1, y1, x2, y2 = self.x[0], self.y[0], self.x[1], self.y[1]
        roi = self.roi
        pts_temp = np.ones(self.pts.shape, dtype=bool)
        d = 0.5 * t
        if x1 == x2: #avoid singularities in slope calculations
            y3, y6 = y1, y1
            y4, y5 = y2, y2
            x3, x4 = x1 - d, x1 - d
            x5, x6 = x1 + d, x1 + d
        elif y1==y2:
            x3, x6 = x1, x1
            x4, x5 = x2, x2
            y3, y4 = y1 - d, y1 - d
            y5, y6 = y1 + d, y1 + d
        else:
            m = (y2 - y1) / (x2 - x1)
            b0 = y2 - m * x2
            b1 = b0 - d * np.sqrt(m ** 2 + 1)
            b2 = b0 + d * np.sqrt(m ** 2 + 1)
            if b2 < b1:
                b1, b2 = b2, b1
            b3 = y1 + x1 / m
            b4 = y2 + x2 / m
            if b4 < b3:
                b3, b4 = b4, b3
            denom = m + (1 / m)
            
            x3 = (b4 - b1) / denom
            x4 = (b3 - b1) / denom
            x5 = (b4 - b2) / denom
            x6 = (b3 - b2) / denom
            
            y3 = m * x3 + b1
            y4 = m * x4 + b1
            y5 = m * x5 + b2
            y6 = m * x6 + b2
            if ((x4 - x5) ** 2 + (y4 - y5) ** 2) > ((x4 - x6) ** 2 + (y4 - y6) **2):
                #correct order of points; keep path from crossing itself
                x5, y5, x6, y6 = x6, y6, x5, y5
                
        verts = np.array([(x3, y3),
                 (x4, y4),
                 (x5, y5),
                 (x6, y6)])
        
        codes = [
        Path.MOVETO,
        Path.LINETO,
        Path.LINETO,
        Path.LINETO]
        
        self.path = Path(verts, codes)
        for i in range(self.pts.shape[0]):
            if roi[0] < i < roi[1]:
                for j in range(self.pts.shape[1]):
                    if roi[2] < j < roi[3]:
                        pts_temp[i, j] = not self.path.contains_point([j, i])
        self.pts *= pts_temp

    def rotate_average(self, fft, mask):
        np.seterr(divide='ignore')
        rotate_avg = np.zeros(fft.shape, dtype='complex64')
        rmax = int(np.ceil(max((mask!=0).sum(0)) / 2))
        r_array = np.zeros((rmax + 1, int(np.ceil(2 * np.pi * rmax) + 10)), dtype='complex64')
        # +1 corrects r_array size to include rmax
        ind = np.zeros(rmax + 1, dtype=int) #indicates how many nonzero elements are in r_array
        r_avg = np.zeros(rmax + 1, dtype='complex64')
        center = ((fft.shape[0] - 1) / 2, (fft.shape[1] - 1) / 2)
        #print(rmax, rotate_avg.shape, r_array.shape, ind.shape, r_avg.shape, center)
        for i in range(fft.shape[0]):
            for j in range(fft.shape[1]):
                if mask[i,j]:
                    r = int(np.sqrt((center[0] - i) ** 2 + (center[1] - j) ** 2))
                    r_array[(r, ind[r])] = fft[i,j]
                    ind[r] += 1 
        r_avg = np.true_divide(r_array.sum(1), ind) #calculate average values for each radius
        for i in range(fft.shape[0]):
            for j in range(fft.shape[1]):
                if not self.pts[i, j] and mask[i, j]:
                    r = int(np.sqrt((center[0] - i) ** 2 + (center[1] - j) ** 2))
                    rotate_avg[i, j] = r_avg[r]           
        return rotate_avg

def filter_demo():                        
    import hyperspy.api as hs
    from tkinter.filedialog import askopenfilename
    
    obj = hs.load(askopenfilename(), signal_type='hologram')
    sb_pos = obj.estimate_sideband_position(sb='lower')
    sb_size = obj.estimate_sideband_size(sb_pos)
    roi = (obj.data.shape[0] // 2 - int(1.2 * sb_size.data[0]), 
           obj.data.shape[0] // 2 + int(1.2 * sb_size.data[0]))
    roi2 = [roi[0], roi[1], roi[0], roi[1]]
    obj_fft = np.fft.fft2(obj.data)
    obj_fft = np.fft.fftshift(obj_fft)
    im = np.log(abs(obj_fft))
    center_b = LineFilter(im, roi2)
    plt.subplot(121)
    plt.imshow(np.multiply(center_b.pts, im))
    plt.subplot(122)
    plt.imshow(np.log(abs(np.fft.ifft2(np.fft.ifftshift(np.multiply(center_b.pts, obj_fft)))) + 1))
    plt.show()    
    
if __name__ == '__main__':
    filter_demo()