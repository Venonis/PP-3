#include <mpi.h>
#include <omp.h>
#include <iostream>

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int threads_count;

    if (world_rank == 0) {
        std::cin >> threads_count;
    }

    MPI_Bcast(&threads_count, 1, MPI_INT, 0, MPI_COMM_WORLD);

#pragma omp parallel num_threads(threads_count)
    {
        int thread_num = omp_get_thread_num();
        int total_threads = threads_count * world_size;

#pragma omp critical
        {
            std::cout << "I am " << thread_num
                << " thread from " << world_rank
                << " process. Number of hybrid threads = " << total_threads
                << std::endl;
        }
    }

    MPI_Finalize();
    return 0;
}
