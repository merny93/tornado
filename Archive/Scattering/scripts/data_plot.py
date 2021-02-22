import numpy as np
import file_tools as ft
import os
import matplotlib.pyplot as plt

def collapse_data(data_list):
    # load the data and collapses the runs in one big array
    from functools import reduce
    return reduce(lambda x,y: x+y["Counts"], data_list, np.zeros_like(data_list[0]["Counts"]))


data_loc = "../data/runs"
run_name = "02_Aluminum"
angles = os.listdir(data_loc)
print("anlges to look at", angles)

for angle in angles:
    data_set = ft.get_data(os.path.join(data_loc, angle, run_name))
    unified_count = collapse_data(data_set)
    plt.plot(unified_count, label=angle)
plt.legend()
plt.savefig("../figures/fast_plot.png")
plt.show()