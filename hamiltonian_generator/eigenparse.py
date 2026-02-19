import datetime

def atomic_number(s):
    assert len(s) < 3
    if s == "H":
        return 1
    elif s == "He":
        return 2
    elif s == "Li":
        return 3
    elif s == "Be":
        return 4
    elif s == "B":
        return 5
    else:
        raise KeyError(f"unknown element \"{s}\"")

def mysortkey(d):
    return (
            d["total"],
            d["Z1"],
            d["Z2"],
            float(d["spacing"]),
            d["basis"],
            d["mapping"],
            )

def main():
    data = list()
    with open("eigen.dat", 'r') as fin:
        t_last = 0.0
        molecule = ""
        spacing = ""
        basis = ""
        active = list(int(0) for x in range(2))
        mapping = ""
        times = dict()
        for line in fin:
            if line[0] != '[':
                continue
            ltokens = line.split(']')
            time_str = ltokens[0][1:]
            message = ltokens[1].strip()
            timestamp = datetime.datetime.strptime(time_str, "%Y-%m-%d %H:%M:%S.%f").timestamp()
            mtokens = message.split(' ')
            if mtokens[0] == "Processing":
                if len(times) > 0:
                    data.append({
                        "Z1": atomic_number(molecule[:molecule.find('-')]),
                        "Z2": atomic_number(molecule[molecule.find('-')+1:]),
                        "molecule": molecule,
                        "spacing": spacing,
                        "basis": basis,
                        "mapping": mapping,
                        "occupied": active[0],
                        "vacant": active[1],
                        "total": sum(active),
                        "2^q": 2**sum(active),
                        "loading": times["loading"],
                        "one-norm": times["one-norm"],
                        "Hamiltonian": times["Hamiltonian"],
                        "largest": times["largest"],
                        "smallest": times["smallest"],
                        })
                t0 = timestamp
                path = message.split('"')[-2]
                file = path.split('/')[-1]
                ftokens = file.split('_')
                molecule = ftokens[0]
                spacing = ftokens[1]
                basis = ftokens[2]
                active = list(int(n) for n in ftokens[3].split('-')[1:])
                mapping = ftokens[4]
                mapping = mapping[:mapping.find('.')]
                new_msg = [
                        f"starting time = {t0}",
                        f"molecule = {molecule}",
                        f"spacing = {spacing}",
                        f"basis = {basis}",
                        f"active_space = {active[0]} + {active[1]} = {sum(active)}",
                        f"mapping = {mapping}",
                        ]
                times.clear()
            else:
                assert mtokens[0] == "--"
                mtokens = mtokens[1:]
                if mtokens[0] == "Beginning":
                    kind = "begin"
                elif mtokens[0] == "Loaded":
                    kind = "loading"
                elif mtokens[0] == "Computed":
                    kind = mtokens[2]
                    if kind[-1] == ':':
                        kind = kind[:-1]
                elif mtokens[0] == "Constructed":
                    kind = "Hamiltonian"
                else:
                    raise ValueError(f"unknown computation: \"{' '.join(mtokens)}\"")
                times[kind] = timestamp - t_last
                t_last = timestamp
    for item in sorted(data, key=mysortkey):
        cells = [
                item["Z1"],
                item["molecule"],
                item["spacing"],
                item["basis"],
                item["mapping"],
                item["occupied"],
                item["vacant"],
                item["total"],
                item["2^q"],
                item["loading"],
                item["one-norm"],
                item["Hamiltonian"],
                item["largest"],
                item["smallest"],
                ]
        print(", ".join(str(column) for column in cells))

if __name__ == "__main__":
    main()
