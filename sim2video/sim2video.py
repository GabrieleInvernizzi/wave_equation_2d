import sys
from pathlib import Path

import imageio as iio
import numpy as np


def main():
    if len(sys.argv) != 2:
        print("Usage: python3 sim2video.py <filename.sim>")
        sys.exit(1)

    sim_name = Path(sys.argv[1]).stem

    print(f"Reading \"{sys.argv[1]}\"...")
    try: 
        f = open(sys.argv[1], "rb")
        sim_data = f.read()
        f.close()
    except IOError:
        print(f"Could not read the file {sys.argv[1]}!")
        sys.exit(2)

    print(f"Rendering \"{sim_name}\"...")

    info, sim_data_raw = sim_data.split(b'\n', maxsplit=1)
    
    dt, rows, cols, steps, fps = info.split(b'-')
    dt, rows, cols, steps, fps = int(dt), int(rows), int(cols), int(steps), int(fps)
    
    if dt == 8:
        sim_data = np.frombuffer(sim_data_raw, np.float64).reshape((steps, rows, cols))
    elif dt == 4:
        sim_data = np.frombuffer(sim_data_raw, np.float32).reshape((steps, rows, cols))
    else:
        print(f"Can't work with double of size = {dt}. Exiting.")
        sys.exit(3)


    min_simdata = np.min(sim_data)
    max_simdata = np.max(sim_data)

    range_simdata = max_simdata - min_simdata;
    if range_simdata < 0.00001:
        print("All data points close to zero. Exiting.")
        sys.exit(0)

    sim_data = np.array(((sim_data - min_simdata) / range_simdata) * 255, dtype=np.uint8)
    
    print(f"Encoding \"{sim_name}\" to avi...")
    w = iio.get_writer(f'{sim_name}.avi', fps=fps)
    for i in range(steps):
        im = iio.core.util.Array(sim_data[i,:,:])
        w.append_data(im)
    w.close()
          
    print(f"Created the video \"{sim_name}.avi\" from {sys.argv[1]}.")



if __name__ == "__main__":
    main()
