#include <stdio.h>

#include "sims.h"

SimConf get_sim_conf() {
    SimConf c = {.filepath = "out.sim",
                 .save_period = 2,
                 .duration = 10,
                 .dt = 0.01,
                 .c = 1.0,
                 .domain_width = 5,
                 .domain_height = 5,
                 .dx = 0.01,
                 .dy = 0.01};

    c.n_steps = (size_t)(c.duration / c.dt);
    c.framerate = 1.0 / (c.dt * c.save_period);
    c.cols = (size_t)(c.domain_width / c.dx);
    c.rows = (size_t)(c.domain_height / c.dy);
    c.tot_cols = c.cols + 2;
    c.tot_rows = c.rows + 2;

    return c;
}

int main() {
    const SimConf c = get_sim_conf();
    // Create the file to save all the frames
    FILE *f = fopen(c.filepath, "wb");
    if (!f) {
        printf("Can't create the file \"%s\". Exiting.\n", c.filepath);
        return 1;
    }

    fprintf(f, "%zu-%zu-%zu-%zu-%u\n", sizeof(double), c.tot_rows, c.tot_cols,
            (size_t)(c.n_steps / c.save_period), c.framerate);

    sim_circ(c, f);

    fclose(f);

    return 0;
}
