#include <mpi.h>
#include <omp.h>
#include <iostream>
#include <iomanip>

int main(int argc, char* argv[]) {
    int rank, size;
    long long num_intervals;
    double step, local_sum = 0.0, pi = 0.0;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0) {
        std::cout << "Set precision (number of intervals): " << std::endl;
        std::cin >> num_intervals;
    }

    MPI_Bcast(&num_intervals, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);

    step = 1.0 / (double)num_intervals;

#pragma omp parallel
    {
        double thread_sum = 0.0;
#pragma omp for schedule(static)
        for (long long i = rank; i < num_intervals; i += size) {
            double x = (i + 0.5) * step;
            thread_sum += 4.0 / (1.0 + x * x);
        }
#pragma omp critical
        {
            local_sum += thread_sum;
        }
    }

    local_sum *= step;

    MPI_Reduce(&local_sum, &pi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        std::cout << "Computed value of pi: " << std::fixed << std::setprecision(15) << pi << std::endl;
    }

    MPI_Finalize();
    return 0;
}
