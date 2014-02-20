#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <stdbool.h>
#include <libgen.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <math.h>


char*
strdup(const char *str)
{
	int n = strlen(str) + 1;
	char *dup = malloc(n);
	if (dup) strcpy(dup, str);
	return dup;
}


void
usage(FILE *f, const char *p)
{
	char *dup = strdup(p);
	char *pname = basename(dup);

	const char *u =
	"usage: %s -w width -h height [-x x] [-y y] in-file out-file.\n"
	"extrapolate optical flow from a lower resolution FLO(W) file to a higher resolution.\n"
	"The output file format (either FLO or FLO) is determined by the input format\n"
	"example: %s -w 512 -h 488 -x 5.0 -y 5.0 small.flow large.flow\n\n"
	"   -H          print this help text\n\n"
	"Extrapolation Control:\n"
	"   -w width    target width\n"
	"   -h height   target height\n"
	"   -x x        extrapolation factor in x-direction\n"
	"   -y y        extrapolation factor in y-direction\n\n"
	"Note that if x and y are not passed, they are determined by input-{w,h}/output-{w,h}\n";

	fprintf(f, u, pname, pname);
	free(dup);
}


bool
test_file(char *f)
{
	struct stat sb;
	if (stat(f, &sb) == -1)
			return false;
	return true;
}

bool
test_files(char *in, char *out)
{
	if (!test_file(in)) {
		fprintf(stderr, "Error: Unavailable file '%s'\n", in);
		return false;
	}
	return true;
}


typedef enum {
	FF_FLO,
	FF_FLOW,
	FF_LAST
} FLOformat;


// header identification strings
static const char *flo_formats[FF_LAST] = {"PIEH", "PIEI"};


typedef struct {
	int w, h;
	FLOformat fmt;
	int nmemb;
	float *data;
} FLOfile;


FLOfile*
read_flow(char *f)
{
	FILE *fd = fopen(f, "r");
	if (!fd) {
		fprintf(stderr, "Error: Could not open file '%s'.\n", f);
		return NULL;
	}

	char *tag = calloc(4, sizeof(char));
	if (fread(tag, sizeof(char), 4, fd) < 4) {
		fprintf(stderr, "Error: Could not read file header of '%s'\n", f);
		fclose(fd);
		return NULL;
	}

	FLOformat fmt = FF_LAST;
	if (!strncmp(tag, flo_formats[FF_FLO], 4)) {
		fmt = FF_FLO;
	}
	else if (!strncmp(tag, flo_formats[FF_FLOW], 4)) {
		fmt = FF_FLOW;
	}
	else {
		fprintf(stderr, "Error: Invalid file type in file '%s'\n", f);
		fclose(fd);
		return NULL;
	}
	free(tag);

	int w;
	if (fread(&w, sizeof(int), 1, fd) < 1) {
		fprintf(stderr, "Error: Could not read width from file '%s'\n", f);
		fclose(fd);
		return NULL;
	}

	int h;
	if (fread(&h, sizeof(int), 1, fd) < 1) {
		fprintf(stderr, "Error: Could not read height from file '%s'\n", f);
		fclose(fd);
		return NULL;
	}

	if (w <= 0 || h <= 0) {
		fprintf(stderr, "Error: Invalid width or height in file '%s'\n", f);
		fclose(fd);
		return NULL;
	}

	int sz = w * h;
	switch (fmt) {
	case FF_FLO:
		sz *= 2;
		break;
	case FF_FLOW:
		sz *= 3;
		break;
	default:
		// should not happen
		break;
	}
	FLOfile *flow = malloc(sizeof(FLOfile));
	flow->nmemb = sz;
	flow->fmt = fmt;
	flow->w = w;
	flow->h = h;
	flow->data = malloc(flow->nmemb * sizeof(float));
	if (fread(flow->data, 1, flow->nmemb * sizeof(float), fd) < flow->nmemb * sizeof(float)) {
		fprintf(stderr, "Error: Incomplete data in file '%s'\n", f);
		fclose(fd);
		free(flow->data);
		free(flow);
		return NULL;
	}

	fclose(fd);
	return flow;
}


int
write_flow(char *f, FLOfile *flo)
{
	FILE *fd = fopen(f, "w");
	if (!fd) {
		fprintf(stderr, "Error: Could not open '%s' for writing\n", f);
		return -1;
	}

	fprintf(fd, "%s", flo_formats[flo->fmt]);
	fwrite(&flo->w, sizeof(int), 1, fd);
	fwrite(&flo->h, sizeof(int), 1, fd);
	fwrite(flo->data, 1, flo->nmemb * sizeof(float), fd);
	fclose(fd);

	return 0;
}


void
free_flow(FLOfile *f)
{
	if (f) free(f->data);
	free(f);
}


void
interpolate2D_weights(float rx, float ry, float w[2][2])
{
	int ix, iy;
	float x, y;

	iy = (int)floor(ry);
	ix = (int)floor(rx);
	y = ry - (float)iy;
	x = rx - (float)ix;

	w[0][0] = (1 - x) * (1 - y);
	w[1][0] = (x) * (1 - y);
	w[0][1] = (1 - x) * (y);
	w[1][1] = (x) * (y);
}


#define _U_IDX(flow,x,y, npb) ((y) * ((flow)->w * (npb)) + (x) * (npb) + 0)
#define _V_IDX(flow,x,y, npb) ((y) * ((flow)->w * (npb)) + (x) * (npb) + 1)
#define _C_IDX(flow,x,y, npb) ((y) * ((flow)->w * (npb)) + (x) * (npb) + 2)
#define _U(flow, x, y, npb) ((flow)->data[_U_IDX((flow),(x),(y),(npb))])
#define _V(flow, x, y, npb) ((flow)->data[_V_IDX((flow),(x),(y),(npb))])
#define _C(flow, x, y, npb) ((flow)->data[_C_IDX((flow),(x),(y),(npb))])


int
extrapolate(const int w, const int h, float sx, float sy, char *in, char *out)
{
	if (!test_files(in, out)) return EXIT_FAILURE;

	FLOfile *in_flo = read_flow(in);
	if (!in_flo) return EXIT_FAILURE;

	FLOfile *out_flo = malloc(sizeof(FLOfile));
	out_flo->fmt = in_flo->fmt;
	int npb = 2; // number of items per block (one block = vx, vy, possibly confidence)
	switch (in_flo->fmt) {
	case FF_FLO:
		npb = 2;
		break;
	case FF_FLOW:
		npb = 3;
		break;
	default:
		// should not happen
		break;
	}
	out_flo->nmemb = npb * w * h;
	out_flo->w = w;
	out_flo->h = h;
	out_flo->data = malloc(out_flo->nmemb * sizeof(float));
	memset(out_flo->data, 0, out_flo->nmemb * sizeof(float));

	// determine exact scalings if none are passed
	if (sx <= 0.0f) sx = (float)w / (float)in_flo->w;
	if (sy <= 0.0f) sy = (float)h / (float)in_flo->h;

	float weight[2][2];
	for (int y = 0; y < h; y++) {
		for (int x = 0; x < w; x++) {
			float rx = (float)x / sx;
			float ry = (float)y / sy;


			int kx1 = (int)floor(rx);
			int ky1 = (int)floor(ry);
			int kx2 = kx1 + 2;
			int ky2 = ky1 + 2;

			interpolate2D_weights(rx, ry, weight);
			int wx, wy, ix, iy;
			float u = 0.0f;
			float v = 0.0f;
			float c = 0.0f;
			int n = 0;
			for (iy = ky1, wy = 0; iy < ky2; iy++, wy++) {
				for (ix = kx1, wx = 0; ix < kx2; ix++, wx++, n++) {
					if (ix >= in_flo->w || iy >= in_flo->h) continue;
					u += weight[wx][wy] * _U(in_flo, ix, iy, npb);
					v += weight[wx][wy] * _V(in_flo, ix, iy, npb);
					if (out_flo->fmt == FF_FLOW)
						c += weight[wx][wy] * _C(in_flo, ix, iy, npb);
				}
			}
			_U(out_flo, x, y, npb) = u;
			_V(out_flo, x, y, npb) = v;
			if (out_flo->fmt == FF_FLOW)
				_C(out_flo, x, y, npb) = c;
		}
	}

	write_flow(out, out_flo);
	free_flow(in_flo);
	free_flow(out_flo);
	return EXIT_SUCCESS;
}


int
main(int argc, char *argv[])
{
	const char *pname = argv[0];
	if (argc <= 1) {
		usage(stderr, pname);
		return EXIT_FAILURE;
	}

	bool wb, hb;
	int w, h, opt;
	float sx, sy;

	wb = hb = false;
	w = h = 0;
	sx = sy = -1.0;
	while ((opt = getopt(argc, argv, "Hw:h:x:y:")) != -1) {
		switch(opt) {
		case 'H':
			usage(stdout, pname);
			return EXIT_SUCCESS;
		case 'h':
			h = atoi(optarg);
			hb = true;
			break;
		case 'w':
			w = atoi(optarg);
			wb = true;
			break;
		case 'x':
			sx = atof(optarg);
			break;
		case 'y':
			sy = atof(optarg);
			break;
		default:
			fprintf(stderr, "Unknown argument '%c'\n", opt);
			usage(stderr, pname);
			return EXIT_FAILURE;
		}
	}

	if (!wb || !hb) {
		fprintf(stderr, "Missing argument for width or height\n");
		return EXIT_FAILURE;
	}

	if (argc - optind != 2) {
		fprintf(stderr, "Expected two arguments (input file, output file) after options.\n");
		return EXIT_FAILURE;
	}

	return extrapolate(w, h, sx, sy, argv[optind], argv[optind+1]);
}
