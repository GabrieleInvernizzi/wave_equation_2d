#include "sim_conf.h"


SimConf get_sim_conf() {
	SimConf c = {
		.filepath = "out.sim",
		.save_period = 2,
		.duration = 10,
		.dt = 0.01,
		.c = 1.0,
		.domain_width = 5,
		.domain_height = 5,
		.dx = 0.01,
		.dy = 0.01
	};

	c.n_steps = (size_t)(c.duration / c.dt);
	c.framerate = 1.0/(c.dt*c.save_period);
	c.cols = (size_t)(c.domain_width / c.dx);
	c.rows = (size_t)(c.domain_height / c.dy);
	c.tot_cols = c.cols + 2;
	c.tot_rows = c.rows + 2;

	return c;
}
