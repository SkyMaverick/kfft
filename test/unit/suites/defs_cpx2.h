#pragma once

#include <kfft.h>

#define FFTX_1011 10
#define FFTY_1011 11
const kfft_cpx fftin_1011[FFTX_1011 * FFTY_1011] = {
    {1.0000, 0.0000},  {2.0000, 0.0000},  {3.0000, 0.0000},  {4.0000, 0.0000},  {5.0000, 0.0000},
    {6.0000, 0.0000},  {7.0000, 0.0000},  {8.0000, 0.0000},  {9.0000, 0.0000},  {10.0000, 0.0000},
    {2.0000, 0.0000},  {3.0000, 0.0000},  {4.0000, 0.0000},  {5.0000, 0.0000},  {6.0000, 0.0000},
    {7.0000, 0.0000},  {8.0000, 0.0000},  {9.0000, 0.0000},  {10.0000, 0.0000}, {11.0000, 0.0000},
    {3.0000, 0.0000},  {4.0000, 0.0000},  {5.0000, 0.0000},  {6.0000, 0.0000},  {7.0000, 0.0000},
    {8.0000, 0.0000},  {9.0000, 0.0000},  {10.0000, 0.0000}, {11.0000, 0.0000}, {12.0000, 0.0000},
    {4.0000, 0.0000},  {5.0000, 0.0000},  {6.0000, 0.0000},  {7.0000, 0.0000},  {8.0000, 0.0000},
    {9.0000, 0.0000},  {10.0000, 0.0000}, {11.0000, 0.0000}, {12.0000, 0.0000}, {13.0000, 0.0000},
    {5.0000, 0.0000},  {6.0000, 0.0000},  {7.0000, 0.0000},  {8.0000, 0.0000},  {9.0000, 0.0000},
    {10.0000, 0.0000}, {11.0000, 0.0000}, {12.0000, 0.0000}, {13.0000, 0.0000}, {14.0000, 0.0000},
    {6.0000, 0.0000},  {7.0000, 0.0000},  {8.0000, 0.0000},  {9.0000, 0.0000},  {10.0000, 0.0000},
    {11.0000, 0.0000}, {12.0000, 0.0000}, {13.0000, 0.0000}, {14.0000, 0.0000}, {15.0000, 0.0000},
    {7.0000, 0.0000},  {8.0000, 0.0000},  {9.0000, 0.0000},  {10.0000, 0.0000}, {11.0000, 0.0000},
    {12.0000, 0.0000}, {13.0000, 0.0000}, {14.0000, 0.0000}, {15.0000, 0.0000}, {16.0000, 0.0000},
    {8.0000, 0.0000},  {9.0000, 0.0000},  {10.0000, 0.0000}, {11.0000, 0.0000}, {12.0000, 0.0000},
    {13.0000, 0.0000}, {14.0000, 0.0000}, {15.0000, 0.0000}, {16.0000, 0.0000}, {17.0000, 0.0000},
    {9.0000, 0.0000},  {10.0000, 0.0000}, {11.0000, 0.0000}, {12.0000, 0.0000}, {13.0000, 0.0000},
    {14.0000, 0.0000}, {15.0000, 0.0000}, {16.0000, 0.0000}, {17.0000, 0.0000}, {18.0000, 0.0000},
    {10.0000, 0.0000}, {11.0000, 0.0000}, {12.0000, 0.0000}, {13.0000, 0.0000}, {14.0000, 0.0000},
    {15.0000, 0.0000}, {16.0000, 0.0000}, {17.0000, 0.0000}, {18.0000, 0.0000}, {19.0000, 0.0000},
    {11.0000, 0.0000}, {12.0000, 0.0000}, {13.0000, 0.0000}, {14.0000, 0.0000}, {15.0000, 0.0000},
    {16.0000, 0.0000}, {17.0000, 0.0000}, {18.0000, 0.0000}, {19.0000, 0.0000}, {20.0000, 0.0000}};

const kfft_cpx fftout_1011[FFTX_1011 * FFTY_1011] = {
    {66.0000, 0.0000},   {77.0000, 0.0000},   {88.0000, 0.0000},   {99.0000, 0.0000},
    {110.0000, 0.0000},  {121.0000, 0.0000},  {132.0000, 0.0000},  {143.0000, 0.0000},
    {154.0000, 0.0000},  {165.0000, 0.0000},  {-5.5000, 18.7313},  {-5.5000, 18.7313},
    {-5.5000, 18.7313},  {-5.5000, 18.7313},  {-5.5000, 18.7313},  {-5.5000, 18.7313},
    {-5.5000, 18.7313},  {-5.5000, 18.7313},  {-5.5000, 18.7313},  {-5.5000, 18.7313},
    {-5.5000, 8.5582},   {-5.5000, 8.5582},   {-5.5000, 8.5582},   {-5.5000, 8.5582},
    {-5.5000, 8.5582},   {-5.5000, 8.5582},   {-5.5000, 8.5582},   {-5.5000, 8.5582},
    {-5.5000, 8.5582},   {-5.5000, 8.5582},   {-5.5000, 4.7658},   {-5.5000, 4.7658},
    {-5.5000, 4.7658},   {-5.5000, 4.7658},   {-5.5000, 4.7658},   {-5.5000, 4.7658},
    {-5.5000, 4.7658},   {-5.5000, 4.7658},   {-5.5000, 4.7658},   {-5.5000, 4.7658},
    {-5.5000, 2.5118},   {-5.5000, 2.5118},   {-5.5000, 2.5118},   {-5.5000, 2.5118},
    {-5.5000, 2.5118},   {-5.5000, 2.5118},   {-5.5000, 2.5118},   {-5.5000, 2.5118},
    {-5.5000, 2.5118},   {-5.5000, 2.5118},   {-5.5000, 0.7908},   {-5.5000, 0.7908},
    {-5.5000, 0.7908},   {-5.5000, 0.7908},   {-5.5000, 0.7908},   {-5.5000, 0.7908},
    {-5.5000, 0.7908},   {-5.5000, 0.7908},   {-5.5000, 0.7908},   {-5.5000, 0.7908},
    {-5.5000, -0.7908},  {-5.5000, -0.7908},  {-5.5000, -0.7908},  {-5.5000, -0.7908},
    {-5.5000, -0.7908},  {-5.5000, -0.7908},  {-5.5000, -0.7908},  {-5.5000, -0.7908},
    {-5.5000, -0.7908},  {-5.5000, -0.7908},  {-5.5000, -2.5118},  {-5.5000, -2.5118},
    {-5.5000, -2.5118},  {-5.5000, -2.5118},  {-5.5000, -2.5118},  {-5.5000, -2.5118},
    {-5.5000, -2.5118},  {-5.5000, -2.5118},  {-5.5000, -2.5118},  {-5.5000, -2.5118},
    {-5.5000, -4.7658},  {-5.5000, -4.7658},  {-5.5000, -4.7658},  {-5.5000, -4.7658},
    {-5.5000, -4.7658},  {-5.5000, -4.7658},  {-5.5000, -4.7658},  {-5.5000, -4.7658},
    {-5.5000, -4.7658},  {-5.5000, -4.7658},  {-5.5000, -8.5582},  {-5.5000, -8.5582},
    {-5.5000, -8.5582},  {-5.5000, -8.5582},  {-5.5000, -8.5582},  {-5.5000, -8.5582},
    {-5.5000, -8.5582},  {-5.5000, -8.5582},  {-5.5000, -8.5582},  {-5.5000, -8.5582},
    {-5.5000, -18.7313}, {-5.5000, -18.7313}, {-5.5000, -18.7313}, {-5.5000, -18.7313},
    {-5.5000, -18.7313}, {-5.5000, -18.7313}, {-5.5000, -18.7313}, {-5.5000, -18.7313},
    {-5.5000, -18.7313}, {-5.5000, -18.7313}};
