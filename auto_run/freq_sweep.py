"""
At fixed position sweet frequeuncy and report amplitude

Set position,
Move to position 

set frequency range and frequency resolution

for each frequency in range
measure amplitude by fitting sine function!
"""

import numpy as np
from time import sleep

class FreqSweep:
    def __init__(self, motors, sound_object):
        """
            Position is tuple r(mm), theta(deg)
            frequency range is array of all frequencies to sweep
            motors is a mcphysics motor object
            sound_object is a mcphysics sound card object
        """
        
        
        self.m = motors
        self.so = sound_object
    
    def move_to_pos(self):
        #uses mcphysics.experiments.drum
        self.m.set_ra(self.pos[0], self.pos[1])
        ##this is a blocking function so i guess its gonna lock the code here....
    
    def do_sweep(self, frequency_range, position, use_soundcard = True):

        self.pos = position
        self.freqs = frequency_range
        #record time
        data_time = 0.1
        #sample rate
        sp = 44000
        #waittime
        ring_down = 0.3
        
        #init list of datas
        data = []

        self.move_to_pos()
        if use_soundcard:
            #init the player
            self.so.button_play.set_checked(True)
            for freq in self.freqs:
                #for each freq
                #turn on the freqeuncy
                self.so.tab_out.settings["Left/Sine"] = freq #this is in hertz
                #wait for ring down
                sleep(ring_down)
                #take the data! this is blocking
                d_run = self.so.api.rec(frames = int(data_time*sp), channels=1, blocking=True)
                data.append(d_run)
            #take the statistical max over the lists
            res =list(map(lambda x: np.mean(x) + 5*np.std(x), data))

            #turn off the player
            self.so.button_play.set_checked(True)
            
            return res
        else:
            print("the adlm doesnt work yet")
            return None
        

