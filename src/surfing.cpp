/*
 * Copyright (C) 2014 Bruno Roy
 *
 * This file is part of Surfing.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * blobserver is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with blobserver.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "surfaceReconstruction.h"

int main(int argc, char** argv)
{
    SurfaceReconstruction surfaceReconstruction(argc, argv);
    surfaceReconstruction.reconstruct();

    return 0;
}

