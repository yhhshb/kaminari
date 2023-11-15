import sys

if __name__ == "__main__":
    if len(sys.argv) != 4: raise RuntimeError("Usage: {} [{}] [{}] [{}]".format(sys.argv[0], "input_filename", "sampling_rate", "output_filename"))
    input_filename = sys.argv[1]
    sampling_rate = int(sys.argv[2])
    output_filename = (sys.argv[3])
    with open(output_filename, "w") as oh:
        with open(input_filename, "r") as ih:
            for line in ih:
                line = line.strip()
                if line and ((hash(line) % sampling_rate) == 0):
                    oh.write(line + "\n")
