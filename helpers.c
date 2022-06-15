#include "helpers.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <fenv.h>

// I'm going to create some data structs to ease the process of calculating some stuff.
const int Gx_kernel[][3] =
    {
        {-1, 0, 1},
        {-2, 0, 2},
        {-1, 0, 1}};

const int Gy_kernel[][3] =
    {
        {-1, -2, -1},
        {0, 0, 0},
        {1, 2, 1}};

// Grid data structure.
typedef struct
{
    int pixels;
    int pixels_w;
    int pixels_h;
    int height_top;
    int height_bottom;
    int width_leftboundary;
    int width_rightboundary;
    RGBTRIPLE coordinates[3][3];
    RGBTRIPLE vector[9];
} __attribute__((__packed__))
GRID;

long averages(int numbers, BYTE vector[numbers]);
void get_grid(GRID *target, int height, int width, RGBTRIPLE image[height][width], int x, int y);
BYTE sobel(long Gx, long Gy);

// Convert image to grayscale
void grayscale(int height, int width, RGBTRIPLE image[height][width])
{
    // Grayscale filter.
    //  Iterate over every triplet changing the value.
    //  For every line of pixels…
    for (int i = 0; i < height; i++)
    {
        // For every triplet pixel in this line.
        for (int j = 0; j < width; j++)
        {
            int average = 0;
            if (image[i][j].rgbtBlue == image[i][j].rgbtGreen && image[i][j].rgbtRed == image[i][j].rgbtBlue)
            {
                average = image[i][j].rgbtBlue;
            }
            else
            {
                // Make the average of all three values.
                // Sum of values in a double.
                int sum = (int)image[i][j].rgbtBlue + (int)image[i][j].rgbtGreen + (int)image[i][j].rgbtRed;

                // Divide by 3 values.
                double result = sum / 3.0;

                // Save the round int
                average = (int)nearbyint(result);
            }
            // Save the new value.
            image[i][j].rgbtBlue = average;
            image[i][j].rgbtGreen = average;
            image[i][j].rgbtRed = average;
        }
    }
    return;
}

// Reflect image horizontally
void reflect(int height, int width, RGBTRIPLE image[height][width])
{
    // Reflect.
    // Logic to verify width size.
    int half;
    if (width % 2 == 0)
    {
        half = width / 2;
    }
    else
    {
        half = rint(width / 2);
    }

    // Iterate.
    for (int i = 0; i < height; i++)
    {
        // Third variable.
        RGBTRIPLE tmpPixel;
        for (int j = 0; j < half; j++)
        {
            // Save current pixel in temp pixel.
            tmpPixel = image[i][j];

            // Move mirror pixel to current location.
            image[i][j] = image[i][width - (j + 1)];

            // Move temp pixel to mirror.
            image[i][width - (j + 1)] = tmpPixel;
        }
    }
    return;
}

// FIXME: Blur image
// I think I fixed it, a for loop was running short on its range. Fixed with <=
// Refactoring blur.
void blur(int height, int width, RGBTRIPLE image[height][width])
{
    RGBTRIPLE blurred[height][width];
    GRID *local = malloc(sizeof(GRID));

    // Blur.
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            get_grid(local, height, width, image, i, j);

            // Arrays of color channels.
            BYTE blue[local->pixels];
            BYTE green[local->pixels];
            BYTE red[local->pixels];

            for (int m = 0; m < local->pixels; m++)
            {
                blue[m] = local->vector[m].rgbtBlue;
                green[m] = local->vector[m].rgbtGreen;
                red[m] = local->vector[m].rgbtRed;
            }

            // Get the average for every color channel.
            BYTE blue_average = (BYTE)averages(local->pixels, blue);
            BYTE green_average = (BYTE)averages(local->pixels, green);
            BYTE red_average = (BYTE)averages(local->pixels, red);

            // Save the average value in the original pixel.
            blurred[i][j].rgbtBlue = blue_average;
            blurred[i][j].rgbtGreen = green_average;
            blurred[i][j].rgbtRed = red_average;
        }
    }

    // Free memory.
    free(local);

    // Write the blurred picture over the original picture.
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            image[i][j] = blurred[i][j];
        }
    }
    return;
}

// Detect edges.
void edges(int height, int width, RGBTRIPLE image[height][width])
{
    // FIXME: Edges.
    // Temporary target image.
    RGBTRIPLE edges[height][width];
    GRID *local = malloc(sizeof(GRID));
    long Gx_values[3][3][3];
    long Gy_values[3][3][3];
    if (local == NULL)
    {
        printf("Memory allocation error.\n");
        return;
    }
    long Gx_sum[3];
    long Gy_sum[3];

    // Iterate over the image.
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {

            // Initialize sums to zero.
            for (int k = 0; k < 3; k++)
            {
                Gx_sum[k] = Gy_sum[k] = 0;
            }

            get_grid(local, height, width, image, i, j);
            for (int gi = 0; gi < 3; gi++)
            {
                for (int gj = 0; gj < 3; gj++)
                {
                    // Multiplications per channel.
                    // First for Gx.
                    Gx_values[0][gi][gj] = local->coordinates[gi][gj].rgbtBlue * Gx_kernel[gi][gj];
                    Gx_values[1][gi][gj] = local->coordinates[gi][gj].rgbtGreen * Gx_kernel[gi][gj];
                    Gx_values[2][gi][gj] = local->coordinates[gi][gj].rgbtRed * Gx_kernel[gi][gj];
                    // Now for Gy.
                    Gy_values[0][gi][gj] = local->coordinates[gi][gj].rgbtBlue * Gy_kernel[gi][gj];
                    Gy_values[1][gi][gj] = local->coordinates[gi][gj].rgbtGreen * Gy_kernel[gi][gj];
                    Gy_values[2][gi][gj] = local->coordinates[gi][gj].rgbtRed * Gy_kernel[gi][gj];

                    // Sums for Gx.
                    Gx_sum[0] += Gx_values[0][gi][gj];
                    Gx_sum[1] += Gx_values[1][gi][gj];
                    Gx_sum[2] += Gx_values[2][gi][gj];
                    // Sums for Gy.
                    Gy_sum[0] += Gy_values[0][gi][gj];
                    Gy_sum[1] += Gy_values[1][gi][gj];
                    Gy_sum[2] += Gy_values[2][gi][gj];
                }
            }
            // Calculate the sobel number per channel.
            edges[i][j].rgbtBlue = sobel(Gx_sum[0], Gy_sum[0]);
            edges[i][j].rgbtGreen = sobel(Gx_sum[1], Gy_sum[1]);
            edges[i][j].rgbtRed = sobel(Gx_sum[2], Gy_sum[2]);
        }
    }
    // Free memory allocated.
    free(local);

    // Finally change the image with the edges map.
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            image[i][j] = edges[i][j];
        }
    }
    return;
}

// Extract the surrounding grid of pixels.
// Out of bounds are 0s.
// Grid contains amount of pixels, boundaries and two arrays of RGBTRIPLE values.
void get_grid(GRID *target, int height, int width, RGBTRIPLE image[height][width], int x, int y)
{
    // Get the value of the surrounding pixels per channel.
    // Find the real ammount of pixels that will exist in the blob.
    if (x == 0 || x == height - 1)
    {
        target->pixels_h = 2;
    }
    else
    {
        target->pixels_h = 3;
    }
    if (y == 0 || y == width - 1)
    {
        target->pixels_w = 2;
    }
    else
    {
        target->pixels_w = 3;
    }

    target->pixels = target->pixels_h * target->pixels_w;

    // Find the boundaries to find the pixels since we cannot index beyond the limits.
    // Pixels outside the limits should be 0 at all channels, but should not be present in the array.
    // These variables are the boundaries coordinates of the pixel box.
    // Height.
    if (x == 0)
    {
        target->height_top = 0;
        target->height_bottom = 1;
    }
    else if (x == height - 1)
    {
        target->height_top = -1;
        target->height_bottom = 0;
    }
    else
    {
        target->height_top = -1;
        target->height_bottom = 1;
    }
    // Width.
    if (y == 0)
    {
        target->width_leftboundary = 0;
        target->width_rightboundary = 1;
    }
    else if (y == width - 1)
    {
        target->width_leftboundary = -1;
        target->width_rightboundary = 0;
    }
    else
    {
        target->width_leftboundary = -1;
        target->width_rightboundary = 1;
    }

    // Make Grid blank.
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            target->coordinates[i][j].rgbtBlue = target->coordinates[i][j].rgbtGreen = target->coordinates[i][j].rgbtRed = (BYTE)0;
        }
    }

    // God help me, this should construct a vector of selected pixels.
    int k = 0;
    do
    {
        int s = target->height_top;
        do
        {
            int t = target->width_leftboundary;
            do
            {
                target->vector[k] = target->coordinates[s + 1][t + 1] = image[x + s][y + t];
                t++;
                k++;
            } while (t <= target->width_rightboundary);
            s++;
        } while (s <= target->height_bottom);
    } while (k < target->pixels);

    return;
}

long averages(int numbers, BYTE vector[numbers])
{
    if (fegetround() != FE_TONEAREST)
    {
        fesetround(FE_TONEAREST);
    }

    double sum = 0;
    for (int i = 0; i < numbers; i++)
    {
        sum += (double)vector[i];
    }
    double result = sum / (double)numbers;
    result = nearbyint(result);
    return (long)result;
}

// Returns the Sobel operation of two given integers as a single color byte.
BYTE sobel(long Gx, long Gy)
{
    // Powers.
    double x = pow(Gx, 2);
    double y = pow(Gy, 2);

    // Square root of the addition.
    double root = sqrt((x + y));
    long result = (long)nearbyint(root);

    switch (result)
    {
        case 255 ... INT64_MAX:
            return (BYTE)255;

        default:
            return (BYTE)result;
    }
}
