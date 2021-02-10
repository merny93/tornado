SOURCE_NAMES = ("Ba-133", "Co-57", "Cs-137", "Na-22")

def csv_generic(fname):
    res = {}
    fields = None
    with open(fname, "r") as csvfile:
        for line in csvfile.readlines():
            if fields is None:
                #init the results
                fields = line.strip().split(",")
                res.update(dict(zip(fields, [[] for _ in range(len(fields))])))
            else:
                vals = [float(x) for x in line.strip().split(",")]
                pairs = dict(zip(fields, vals))
                for key in pairs:
                    res[key].append(pairs[key])   
    return res

def _reader(fname):
    print(fname)
    if fname.split(".")[-1]=="Chn":
        res = read_chn(fname)
    else:
        res = read_csv(fname)
    return res

def read_chn(fname):
    from chn2csv import load_chn as _load_chn
    import numpy as np
    res = _load_chn(fname, delimiter=",")
    res_format = {}
    res_format["Channel"] = np.array(res["Channel"], dtype= int)
    res_format["Counts"] = np.array(res["Counts"])
    for key in res.hkeys:
        res_format[key] = res.h(key)
    return res_format


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

def get_run_names(folder_loc, source_name=None):
    """
        returns the files with path that fit the name (delimited by _)
    """
    import os
    all_files = os.listdir(folder_loc)
    if source_name is None:
        return [os.path.join(folder_loc, f) for f in all_files]
    paths_fit = [os.path.join(folder_loc, f) for f in all_files if f.split("_")[0] == source_name]
    return paths_fit

def get_data(loc, source_name=None):
    to_search = get_run_names(loc, source_name=source_name)
    data_set = [_reader(f) for f in to_search]
    return data_set

def read_lit(loc):
    keys = ["energy", "energy_pm", "intensity", "intensity_pm", "dose", "dose_pm"]
    with open(loc, "r") as lit:
        res = {}
        for line in lit.readlines():
            line = line.strip()
            if line[0] =="#":
                continue
            if len(line.split()) == 1:
                current = line
                res[current] = []
                continue
            res[current].append(dict(zip(keys, [float(x) for x in line.split() if x != "%"])))
    return res


# if __name__=="__main__":
#     d = read_csv("Ba-133_Calibration_000.csv")
#     print(type(d["Channel"]))

# if __name__=="__main__":
#     d = read_lit("../calibration_data/litterature")
#     print(d)

if __name__=="__main__":
    d = read_chn("test.Chn")
    print(d.keys())