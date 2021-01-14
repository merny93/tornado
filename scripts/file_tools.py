
def read_csv(fname):
    assert(fname.split(".")[-1] == "csv")
    with open(fname, "r") as csvfile:
        res = {}
        data_ar = [[],[]]
        header_done = False
        for line in csvfile.readlines():
            if line == "\n":
                header_done = True
                continue
            elif header_done:
                for idx,v in enumerate(line.split(",")):
                    data_ar[idx].append(int(v))
            else:
                #doing header
                res[line.split(",")[0]] = line.split(",")[-1]

