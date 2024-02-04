#include <stdio.h>

#include "sim_conf.h"
#include "log.h"
#include "timings.h"
#include "sim.h"


int main(int argc, char** argv) {
    const SimConf c = get_sim_conf(argc, argv);
    // Create the file to save all the frames
    FILE *f = fopen(c.filepath, "wb");
    if (!f) {
        LOGF("Can't create the file \"%s\". Exiting.", c.filepath);
        return 1;
    }

    fprintf(f, "%zu-%zu-%zu-%zu-%u\n", sizeof(double), c.tot_n_cells, c.tot_n_cells,
            (size_t)(c.n_steps / c.save_period), c.framerate);

    sim(c, f);

    fclose(f);

    SAVE_TIMINGS;

    return 0;
}
