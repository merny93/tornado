
def read_csv(fname):
    """
        read in the csv and output dictionary with all meta data
        actual data is in numpy arrays
    """
    import numpy as np
    assert(fname.split(".")[-1] == "csv")
    with open(fname, "r") as csvfile:
        res = {}
        header_done = False
        for line in csvfile.readlines():
            line = line.strip()
            if line == "":
                header_done = True
                continue
            elif header_done:
                try:
                    val_list = [int(v) for v in line.split(",")]
                    for n, f in enumerate(fs):
                        res[f].append(val_list[n])
                except:
                    fs = line.split(",")
                    for f in fs:
                        res[f] = []

            else:
                #doing header
                res[line.split(",")[0]] = line.split(",")[-1]
    for f in fs:
        res[f] = np.array(res[f], dtype=int)
    return res

if __name__=="__main__":
    d = read_csv("Ba-133_Calibration_000.csv")
    print(type(d["Channel"]))