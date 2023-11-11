#include "helpers.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <pthread.h>

#define CONTOUR_CONFIG_COUNT    16
#define FILENAME_MAX_SIZE       50
#define STEP                    8
#define SIGMA                   200
#define RESCALE_X               2048
#define RESCALE_Y               2048

#define CLAMP(v, min, max) if(v < min) { v = min; } else if(v > max) { v = max; }
#define min(a, b) ((a) < (b) ? (a) : (b))

// Creates a map between the binary configuration (e.g. 0110_2) and the corresponding pixels
// that need to be set on the output image. An array is used for this map since the keys are
// binary numbers in 0-15. Contour images are located in the './contours' directory.
typedef struct {
    ppm_image *image;
    ppm_image *new_image;

    int nr_th;
    int thread_id;
    pthread_barrier_t *barrier;

    int step_x;
    int step_y;

    unsigned char **grid;
    ppm_image **contour_map;
} parametrii;




ppm_image **init_contour_map() {
    ppm_image **map = (ppm_image **)malloc(CONTOUR_CONFIG_COUNT * sizeof(ppm_image *));
    if (!map) {
        fprintf(stderr, "Unable to allocate memory\n");
        exit(1);
    }

    for (int i = 0; i < CONTOUR_CONFIG_COUNT; i++) {
        char filename[FILENAME_MAX_SIZE];
        sprintf(filename, "./contours/%d.ppm", i);
        map[i] = read_ppm(filename);
    }

    return map;
}

// Updates a particular section of an image with the corresponding contour pixels.
// Used to create the complete contour image.
void update_image(ppm_image *image, ppm_image *contour, int x, int y) {
    for (int i = 0; i < contour->x; i++) {
        for (int j = 0; j < contour->y; j++) {
            int contour_pixel_index = contour->x * i + j;
            int image_pixel_index = (x + i) * image->y + y + j;

            image->data[image_pixel_index].red = contour->data[contour_pixel_index].red;
            image->data[image_pixel_index].green = contour->data[contour_pixel_index].green;
            image->data[image_pixel_index].blue = contour->data[contour_pixel_index].blue;
        }
    }
}

// Corresponds to step 1 of the marching squares algorithm, which focuses on sampling the image.
// Builds a p x q grid of points with values which can be either 0 or 1, depending on how the
// pixel values compare to the `sigma` reference value. The points are taken at equal distances
// in the original image, based on the `step_x` and `step_y` arguments.
unsigned char **sample_grid(ppm_image *image, int step_x, int step_y, unsigned char sigma) {
    int p = min(image->x, RESCALE_X) / step_x;
    int q = min(image->y, RESCALE_Y) / step_y;

    unsigned char **grid = (unsigned char **)malloc((p + 1) * sizeof(unsigned char*));
    if (!grid) {
        fprintf(stderr, "Unable to allocate memory\n");
        exit(1);
    }

    for (int i = 0; i <= p; i++) {
        grid[i] = (unsigned char *)malloc((q + 1) * sizeof(unsigned char));
        if (!grid[i]) {
            fprintf(stderr, "Unable to allocate memory\n");
            exit(1);
        }
    }

    // for (int i = 0; i < p; i++) {
    //     for (int j = 0; j < q; j++) {
    //         ppm_pixel curr_pixel = image->data[i * step_x * image->y + j * step_y];

    //         unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

    //         if (curr_color > sigma) {
    //             grid[i][j] = 0;
    //         } else {
    //             grid[i][j] = 1;
    //         }
    //     }
    // }
    // grid[p][q] = 0;

    // // last sample points have no neighbors below / to the right, so we use pixels on the
    // // last row / column of the input image for them
    // for (int i = 0; i < p; i++) {
    //     ppm_pixel curr_pixel = image->data[i * step_x * image->y + image->x - 1];

    //     unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

    //     if (curr_color > sigma) {
    //         grid[i][q] = 0;
    //     } else {
    //         grid[i][q] = 1;
    //     }
    // }
    // for (int j = 0; j < q; j++) {
    //     ppm_pixel curr_pixel = image->data[(image->x - 1) * image->y + j * step_y];

    //     unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

    //     if (curr_color > sigma) {
    //         grid[p][j] = 0;
    //     } else {
    //         grid[p][j] = 1;
    //     }
    // }

    return grid;
}

// Corresponds to step 2 of the marching squares algorithm, which focuses on identifying the
// type of contour which corresponds to each subgrid. It determines the binary value of each
// sample fragment of the original image and replaces the pixels in the original image with
// the pixels of the corresponding contour image accordingly.
void march(int thread_id, int nr_th, ppm_image *image, unsigned char **grid, ppm_image **contour_map, int step_x, int step_y) {
    int p = image->x / step_x;
    int q = image->y / step_y;

    int start = thread_id * (double)p / nr_th;
    int end = min((thread_id + 1) * (double)p / nr_th, p);

    for (int i = start; i < end; i++) {
        for (int j = 0; j < q; j++) {
            unsigned char k = 8 * grid[i][j] + 4 * grid[i][j + 1] + 2 * grid[i + 1][j + 1] + 1 * grid[i + 1][j];
            update_image(image, contour_map[k], i * step_x, j * step_y);
        }
    }
}

// Calls `free` method on the utilized resources.
void free_resources(ppm_image *image, ppm_image **contour_map, unsigned char **grid, int step_x) {
    for (int i = 0; i < CONTOUR_CONFIG_COUNT; i++) {
        free(contour_map[i]->data);
        free(contour_map[i]);
    }
    free(contour_map);

    for (int i = 0; i <= image->x / step_x; i++) {
        free(grid[i]);
    }
    free(grid);

    free(image->data);
    free(image);
}

ppm_image *rescale_image(ppm_image *image) {
    // uint8_t sample[3];

    // we only rescale downwards
    if (image->x <= RESCALE_X && image->y <= RESCALE_Y) {
        return image;
    }

    // alloc memory for image
    ppm_image *new_image = (ppm_image *)malloc(sizeof(ppm_image));
    if (!new_image) {
        fprintf(stderr, "Unable to allocate memory\n");
        exit(1);
    }
    new_image->x = RESCALE_X;
    new_image->y = RESCALE_Y;

    new_image->data = (ppm_pixel*)malloc(new_image->x * new_image->y * sizeof(ppm_pixel));
    if (!new_image) {
        fprintf(stderr, "Unable to allocate memory\n");
        exit(1);
    }


    // use bicubic interpolation for scaling
    // for (int i = 0; i < new_image->x; i++) {
    //     for (int j = 0; j < new_image->y; j++) {
    //         float u = (float)i / (float)(new_image->x - 1);
    //         float v = (float)j / (float)(new_image->y - 1);
    //         sample_bicubic(image, u, v, sample);

    //         new_image->data[i * new_image->y + j].red = sample[0];
    //         new_image->data[i * new_image->y + j].green = sample[1];
    //         new_image->data[i * new_image->y + j].blue = sample[2];
    //     }
    // }

    // free(image->data);
    // free(image);

    return new_image;
}
void *f(void *arg)
{
	parametrii *par = (parametrii *)arg;

    ppm_image *image = par->image;
    ppm_image *new_image = par->new_image;

    int nr_th = par->nr_th;
    int thread_id = par->thread_id;
    pthread_barrier_t *barrier = par->barrier;

    int step_x = par->step_x;
    int step_y = par->step_y;
    unsigned char sigma = SIGMA;

    unsigned char **grid = par->grid;
    ppm_image **contour_map = par->contour_map;
    //int ok=0;
    ppm_image *picture;

    if(image->x > RESCALE_X || image->y > RESCALE_Y)
    {
        //ok=1;
        uint8_t sample[3];
        int start = thread_id * (double)new_image->x / nr_th;
        int end = min((thread_id + 1) * (double)new_image->x / nr_th, new_image->x);

        // use bicubic interpolation for scaling
    for (int i = start; i < end; i++) {
        for (int j = 0; j < new_image->y; j++) {
            float u = (float)i / (float)(new_image->x - 1);
            float v = (float)j / (float)(new_image->y - 1);
            sample_bicubic(image, u, v, sample);

            new_image->data[i * new_image->y + j].red = sample[0];
            new_image->data[i * new_image->y + j].green = sample[1];
            new_image->data[i * new_image->y + j].blue = sample[2];
        }
     }
    picture = new_image;

    }else{
        picture = image;
    }
    
    pthread_barrier_wait(barrier);

    
    int p = picture->x / step_x;
    int q = picture->y / step_y;

    int start = thread_id * (double)p / nr_th;
    int end = min((thread_id + 1) * (double)p / nr_th, p);

    for (int i = start; i < end; i++) {
        for (int j = 0; j < q; j++) {
            ppm_pixel curr_pixel = picture->data[i * step_x * picture->y + j * step_y];

            unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

            if (curr_color > sigma) {
                grid[i][j] = 0;
            } else {
                grid[i][j] = 1;
            }
        }
    }
    grid[p][q] = 0;

    // last sample points have no neighbors below / to the right, so we use pixels on the
    // last row / column of the input image for them
    for (int i = start; i < end; i++) {
        ppm_pixel curr_pixel = picture->data[i * step_x * picture->y + picture->x - 1];

        unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

        if (curr_color > sigma) {
            grid[i][q] = 0;
        } else {
            grid[i][q] = 1;
        }
    }

    start = thread_id * (double)q / nr_th;
    end = min((thread_id + 1) * (double)q / nr_th, q);
    for (int j = start; j < end; j++) {
        ppm_pixel curr_pixel = picture->data[(picture->x - 1) * picture->y + j * step_y];

        unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

        if (curr_color > sigma) {
            grid[p][j] = 0;
        } else {
            grid[p][j] = 1;
        }
    }

    pthread_barrier_wait(barrier);
    // 3. March the squares
    march(thread_id, nr_th ,picture, grid, contour_map, step_x, step_y);

	pthread_exit(NULL);
}

int main(int argc, char *argv[]) {
    if (argc < 4) {
        fprintf(stderr, "Usage: ./tema1 <in_file> <out_file> <P>\n");
        return 1;
    }
    int ok =0;

    ppm_image *image = read_ppm(argv[1]);

    int step_x = STEP;
    int step_y = STEP;

    // 0. Initialize contour map
    ppm_image **contour_map = init_contour_map();

    ppm_image *picture;
    // 1. Rescale the image
    ppm_image *scaled_image = rescale_image(image);
    if(image->x > RESCALE_X || image->y > RESCALE_Y) {
        picture = scaled_image;

    }else{
        picture = image;
    }

    // 2. Sample the grid
    unsigned char **grid = sample_grid(picture, step_x, step_y, SIGMA);

    int nr_th = atoi(argv[3]);
    pthread_t threads[nr_th];
    
    int r, i;

    pthread_barrier_t barrier;
    pthread_barrier_init(&barrier, NULL, nr_th);

    parametrii para[nr_th];
    for (i = 0; i < nr_th; i++) {
        para[i].image = image;
        para[i].new_image = scaled_image;

        para[i].nr_th = nr_th;
        para[i].thread_id = i;		

        para[i].barrier = &barrier;
        para[i].step_x = step_x;
        para[i].step_y = step_y;
        para[i].grid = grid;
        para[i].contour_map = contour_map;

	}

    for (i = 0; i < nr_th; i++) {
		r = pthread_create(&threads[i], NULL, f, &para[i]);

		if (r) {
			printf("Eroare la crearea thread-ului %d\n", i);
			exit(-1);
		}
	}
    for (i = 0; i < nr_th; i++) {
		r = pthread_join(threads[i], NULL);

		if (r) {
			printf("Eroare la asteptarea thread-ului %d\n", i);
			exit(-1);
		}
	}

    // 4. Write output
    write_ppm(picture, argv[2]);

    if(picture == image) {
        free(scaled_image->data);
        free(scaled_image);
    }else{
        free(image->data);
        free(image);
    }
    free_resources(picture, contour_map, grid, step_x);

    
    pthread_barrier_destroy(&barrier);

    return 0;
}