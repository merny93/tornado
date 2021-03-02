"""
At fixed position sweet frequeuncy and report amplitude

Set position,
Move to position 

set frequency range and frequency resolution

for each frequency in range
measure amplitude by fitting sine function!
"""

import numpy as np
from time import sleep, time
import matplotlib.pyplot as plt


class FreqSweep:
    def __init__(self, position, frequency_range, motors, sound_object):
        """
            Position is tuple r(mm), theta(deg)
            frequency range is array of all frequencies to sweep
            motors is a mcphysics motor object
            sound_object is a mcphysics sound card object
        """
        
        self.pos = position
        self.freqs = frequency_range
        self.m = motors
        self.so = sound_object
    
    def move_to_pos(self):
        #uses mcphysics.experiments.drum
        self.m.set_ra(self.pos[0], self.pos[1])
        ##this is a blocking function so i guess its gonna lock the code here....
    
    def do_sweep(self, use_soundcard = True):
        #record time
        data_time = 0.2
        #sample rate
        sp = 44000
        #waittime
        ring_down = 1
        
        #init list of datas
        data = []

        self.move_to_pos()
        if use_soundcard:
            for freq in self.freqs:
                print("working on freqency", freq)
                
                #for each freq
                #turn on the freqeuncy
                self.so.tab_out.settings["Left/Sine"] = freq #this is in hertz
                #init the player
                self.so.button_play.set_checked(True)
                
                #wait for ring down
                self.so.window.sleep(ring_down)
                
                #Set itterations to 1
                self.so.tab_in.settings['Iterations'] = 1
                
                #click the record button
                self.so.button_record(True)
                
                #wait for it to finish
                #but because the gui is threaded is has weird race condtions
                #add a timeout
                t = time()
                while self.so.button_record() and time()-t < 0.5: 
                    self.so.window.sleep()
                #if the timout happend
                if self.so.button_record():
                    print("timout on the gui")
                    self.so.button_record(False)
                #grab the data
                d_run = self.so.tab_in.plot_raw["Left"]
                #print(d_run)
                data.append(d_run)
                
                #turn off the player?
                self.so.button_play.set_checked(False)
                
                
            #take the statistical max over the lists
            res =list(map(lambda x: np.max(x), data))

            #turn off the player
            self.so.button_play.set_checked(False)
            
            return res
        else:
            print("the adlm doesnt work yet")
            return None


