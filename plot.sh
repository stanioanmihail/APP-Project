cd MPI/
gnuplot  -persist -e "set xlabel 'Number of Threads'; set ylabel 'Time(seconds)'; set title 'MPI'; set term png; set output 'mpi.png'; plot 'alps_x16.mpirez' using 1:2 title 'Alps x16' with linespoints, 'alps_x18.mpirez' using 1:2 title 'Alps x18' with linespoints, 'alps_x20.mpirez' using 1:2 title 'Alps x20' with linespoints"
cd ../

#cd omp/
#gnuplot  -persist -e "set xlabel 'Number of Threads'; set ylabel 'Time(seconds)'; set title 'OpenMP'; set term png; set output 'omp.png'; plot 'alps_x16.rez' using 1:2 title 'Alps x16' with linespoints, 'alps_x18.rez' using 1:2 title 'Alps x18' with linespoints, 'alps_x20' using 1:2 title 'Alps x20' with linespoints"
#cd ../

#cd PThreads/
#gnuplot  -persist -e "set xlabel 'Number of Threads'; set ylabel 'Time(seconds)'; set title 'PThreads'; set term png; set output 'pthreads.png'; plot 'alps_x16.rez' using 1:2 title 'Alps x16' with linespoints, 'alps_x18.rez' using 1:2 title 'Alps x18' with linespoints, 'alps_x20' using 1:2 title 'Alps x20' with linespoints"
#cd ../

#cd omp_vs_mpi/
#gnuplot  -persist -e "set xlabel 'Number of Threads'; set ylabel 'Time(seconds)'; set title 'OpenMP vs MPI'; set term png; set output 'omp_mpi.png'; plot 'omp500000_3.rez' using 1:2 title '[OMP] 500,000 points' with linespoints, 'omp1000000_3.rez' using 1:2 title '[OMP] 1,000,000 points' with linespoints, 'omp3000000_3.rez' using 1:2 title '[OMP] 3,000,000 points' with linespoints, 'omp5000000_3.rez' using 1:2 title '[OMP] 5,000,000 points' with linespoints, 'mpi500000_3.rez' using 1:2 title '[MPI] 500,000 points' with linespoints, 'mpi1000000_3.rez' using 1:2 title '[MPI] 1,000,000 points' with linespoints, 'mpi3000000_3.rez' using 1:2 title '[MPI] 3,000,000 points' with linespoints, 'mpi5000000_3.rez' using 1:2 title '[MPI] 5,000,000 points' with linespoints"
#cd ../

#cd omp_vs_pthreads/
#gnuplot  -persist -e "set xlabel 'Number of Threads'; set ylabel 'Time(seconds)'; set title 'OpenMP vs PThreads'; set term png; set output 'omp_pthreads.png'; plot 'omp500000_3.rez' using 1:2 title '[OMP] 500,000 points' with linespoints, 'omp1000000_3.rez' using 1:2 title '[OMP] 1,000,000 points' with linespoints, 'omp3000000_3.rez' using 1:2 title '[OMP] 3,000,000 points' with linespoints, 'omp5000000_3.rez' using 1:2 title '[OMP] 5,000,000 points' with linespoints, 'pthreads500000_3.rez' using 1:2 title '[PThreads] 500,000 points' with linespoints, 'pthreads1000000_3.rez' using 1:2 title '[PThreads] 1,000,000 points' with linespoints, 'pthreads3000000_3.rez' using 1:2 title '[PThreads] 3,000,000 points' with linespoints, 'pthreads5000000_3.rez' using 1:2 title '[PThreads] 5,000,000 points' with linespoints"
#cd ../

#cd mpi_vs_pthreads/
#gnuplot  -persist -e "set xlabel 'Number of Threads'; set ylabel 'Time(seconds)'; set title 'MPI vs PThreads'; set term png; set output 'mpi_pthreads.png'; plot 'mpi500000_3.rez' using 1:2 title '[MPI] 500,000 points' with linespoints, 'mpi1000000_3.rez' using 1:2 title '[MPI] 1,000,000 points' with linespoints, 'mpi3000000_3.rez' using 1:2 title '[MPI] 3,000,000 points' with linespoints, 'mpi5000000_3.rez' using 1:2 title '[MPI] 5,000,000 points' with linespoints, 'pthreads500000_3.rez' using 1:2 title '[PThreads] 500,000 points' with linespoints, 'pthreads1000000_3.rez' using 1:2 title '[PThreads] 1,000,000 points' with linespoints, 'pthreads3000000_3.rez' using 1:2 title '[PThreads] 3,000,000 points' with linespoints, 'pthreads5000000_3.rez' using 1:2 title '[PThreads] 5,000,000 points' with linespoints"
#cd ../

#cd omp_mpi_pthreads/
#gnuplot  -persist -e "set xlabel 'Number of Threads'; set ylabel 'Time(seconds)'; set title 'MPI vs PThreads vs OpenMP'; set term png; set output 'mpi_pthreads.png'; plot 'mpi500000_3.rez' using 1:2 title '[MPI] 500,000 points' with linespoints, 'mpi1000000_3.rez' using 1:2 title '[MPI] 1,000,000 points' with linespoints, 'mpi3000000_3.rez' using 1:2 title '[MPI] 3,000,000 points' with linespoints, 'mpi5000000_3.rez' using 1:2 title '[MPI] 5,000,000 points' with linespoints, 'pthreads500000_3.rez' using 1:2 title '[PThreads] 500,000 points' with linespoints, 'pthreads1000000_3.rez' using 1:2 title '[PThreads] 1,000,000 points' with linespoints, 'pthreads3000000_3.rez' using 1:2 title '[PThreads] 3,000,000 points' with linespoints, 'pthreads5000000_3.rez' using 1:2 title '[PThreads] 5,000,000 points' with linespoints, 'omp500000_3.rez' using 1:2 title '[OMP] 500,000 points' with linespoints, 'omp1000000_3.rez' using 1:2 title '[OMP] 1,000,000 points' with linespoints, 'omp3000000_3.rez' using 1:2 title '[OMP] 3,000,000 points' with linespoints, 'omp5000000_3.rez' using 1:2 title '[OMP] 5,000,000 points' with linespoints"
#cd ../
