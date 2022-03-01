#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <pthread.h>
#include "genetic_algorithm.h"

typedef struct thread_argument {
	pthread_barrier_t *barrier;
    individual *current_generation;
	individual *next_generation;
	const sack_object *objects;
	int object_count;
	int generations_count;
	int sack_capacity;
	int thread_id;
	int nr_threads;
} thread_argument;

typedef struct thread_arg_with_id {
	thread_argument *ta;
	int thread_id;
} thread_arg_with_id;


int min(int num1, int num2) 
{
    return (num1 > num2 ) ? num2 : num1;
}

int read_input(sack_object **objects, int *object_count, int *sack_capacity, int *generations_count, int *nr_threads, int argc, char *argv[])
{
	FILE *fp;

	if (argc < 4) {
		fprintf(stderr, "Usage:\n\t./tema1 in_file generations_count\n");
		return 0;
	}

	fp = fopen(argv[1], "r");
	if (fp == NULL) {
		return 0;
	}

	if (fscanf(fp, "%d %d", object_count, sack_capacity) < 2) {
		fclose(fp);
		return 0;
	}

	if (*object_count % 10) {
		fclose(fp);
		return 0;
	}

	sack_object *tmp_objects = (sack_object *) calloc(*object_count, sizeof(sack_object));

	for (int i = 0; i < *object_count; ++i) {
		if (fscanf(fp, "%d %d", &tmp_objects[i].profit, &tmp_objects[i].weight) < 2) {
			free(objects);
			fclose(fp);
			return 0;
		}
	}

	fclose(fp);

	*generations_count = (int) strtol(argv[2], NULL, 10);
	
	if (*generations_count == 0) {
		free(tmp_objects);

		return 0;
	}

	*nr_threads = (int) strtol(argv[3], NULL, 10);

	*objects = tmp_objects;

	return 1;
}

void print_objects(const sack_object *objects, int object_count)
{
	for (int i = 0; i < object_count; ++i) {
		printf("%d %d\n", objects[i].weight, objects[i].profit);
	}
}

void print_generation(const individual *generation, int limit)
{
	for (int i = 0; i < limit; ++i) {
		for (int j = 0; j < generation[i].chromosome_length; ++j) {
			printf("%d ", generation[i].chromosomes[j]);
		}

		printf("\n%d - %d\n", i, generation[i].fitness);
	}
}

void print_best_fitness(const individual *generation)
{
	printf("%d\n", generation[0].fitness);
}

void compute_fitness_function(const sack_object *objects, individual *generation, int object_count, int sack_capacity)
{
	int weight;
	int profit;

	for (int i = 0; i < object_count; ++i) {
		weight = 0;
		profit = 0;

		for (int j = 0; j < generation[i].chromosome_length; ++j) {
			if (generation[i].chromosomes[j]) {
				weight += objects[j].weight;
				profit += objects[j].profit;
			}
		}

		generation[i].fitness = (weight <= sack_capacity) ? profit : 0;
	}
}

void compute_fitness_function_threaded(const sack_object *objects, individual *generation, int object_count, int sack_capacity, int ID, int P)
{
	int weight;
	int profit;
	int start = ID * (double)object_count / P;
	int end = min((ID + 1) * (double)object_count / P, object_count);

	for (int i = start; i < end; ++i) {
		weight = 0;
		profit = 0;

		for (int j = 0; j < generation[i].chromosome_length; ++j) {
			if (generation[i].chromosomes[j]) {
				weight += objects[j].weight;
				profit += objects[j].profit;
			}
		}

		generation[i].fitness = (weight <= sack_capacity) ? profit : 0;
	}
}

int cmpfunc(const void *a, const void *b)
{
	individual *first = (individual *) a;
	individual *second = (individual *) b;

	int res = second->fitness - first->fitness; // decreasing by fitness
	if (res == 0) {

		res = first->nr_chromosomes_1[0] - second->nr_chromosomes_1[0]; // increasing by number of objects in the sack
		
		if (res == 0) {
			return second->index - first->index;
		}
	}

	return res;
}

void mutate_bit_string_1(const individual *ind, int generation_index)
{
	int i, mutation_size;
	int step = 1 + generation_index % (ind->chromosome_length - 2);

	if (ind->index % 2 == 0) {
		// for even-indexed individuals, mutate the first 40% chromosomes by a given step
		mutation_size = ind->chromosome_length * 4 / 10;
		for (i = 0; i < mutation_size; i += step) {
			ind->chromosomes[i] = 1 - ind->chromosomes[i];
			ind->nr_chromosomes_1[0] += ind->chromosomes[i] * 2 - 1;
		}
	} else {
		// for even-indexed individuals, mutate the last 80% chromosomes by a given step
		mutation_size = ind->chromosome_length * 8 / 10;
		for (i = ind->chromosome_length - mutation_size; i < ind->chromosome_length; i += step) {
			ind->chromosomes[i] = 1 - ind->chromosomes[i];
			ind->nr_chromosomes_1[0] += ind->chromosomes[i] * 2 - 1;
		}
	}
}

void mutate_bit_string_2(const individual *ind, int generation_index)
{
	int step = 1 + generation_index % (ind->chromosome_length - 2);

	// mutate all chromosomes by a given step
	for (int i = 0; i < ind->chromosome_length; i += step) {
		ind->chromosomes[i] = 1 - ind->chromosomes[i];
		ind->nr_chromosomes_1[0] += ind->chromosomes[i] * 2 - 1;
	}
}

void crossover(individual *parent1, individual *child1, int generation_index)
{
	individual *parent2 = parent1 + 1;
	individual *child2 = child1 + 1;
	int count = 1 + generation_index % parent1->chromosome_length;

	memcpy(child1->chromosomes, parent1->chromosomes, count * sizeof(int));
	memcpy(child1->chromosomes + count, parent2->chromosomes + count, (parent1->chromosome_length - count) * sizeof(int));

	memcpy(child2->chromosomes, parent2->chromosomes, count * sizeof(int));
	memcpy(child2->chromosomes + count, parent1->chromosomes + count, (parent1->chromosome_length - count) * sizeof(int));

	for(int i = 0; i < child1->chromosome_length; i++) {
		child1->nr_chromosomes_1[0] += child1->chromosomes[i];
	}
	for(int i = 0; i < child2->chromosome_length; i++) {
		child2->nr_chromosomes_1[0] += child2->chromosomes[i];
	}
}

void copy_individual(const individual *from, const individual *to)
{
	memcpy(to->chromosomes, from->chromosomes, from->chromosome_length * sizeof(int));
	to->nr_chromosomes_1[0] = from->nr_chromosomes_1[0];
}

void free_generation(individual *generation)
{
	int i;

	for (i = 0; i < generation->chromosome_length; ++i) {
		free(generation[i].chromosomes);
		free(generation[i].nr_chromosomes_1);
		generation[i].chromosomes = NULL;
		generation[i].fitness = 0;
	}
}

void free_generation_threaded(individual *generation, int ID, int P)
{
	int i;
	int start = ID * (double)generation->chromosome_length / P;
	int end = min((ID + 1) * (double)generation->chromosome_length / P, generation->chromosome_length);

	for (i = start; i < end; ++i) {
		free(generation[i].chromosomes);
		generation[i].chromosomes = NULL;
		generation[i].fitness = 0;
	}
}

void *f(void *arg) {
	thread_arg_with_id a =*(thread_arg_with_id *) arg;
	individual *tmp = NULL;

	int N = a.ta->object_count;
	int P = a.ta->nr_threads;
	int ID = a.thread_id;
	int cursor, count, start, end;

	start = ID * (double)N / P;
	end = min((ID + 1) * (double)N / P, N);

	// set initial generation (composed of object_count individuals with a single item in the sack)
	for (int i = start; i < end; ++i) {
		a.ta->current_generation[i].fitness = 0;
		a.ta->current_generation[i].chromosomes = (int*) calloc(N, sizeof(int));
		a.ta->current_generation[i].chromosomes[i] = 1;
		a.ta->current_generation[i].index = i;
		a.ta->current_generation[i].chromosome_length = N;
		a.ta->current_generation[i].nr_chromosomes_1 = (int*) calloc(1, sizeof(int));
		a.ta->current_generation[i].nr_chromosomes_1[0] = 1;

		a.ta->next_generation[i].fitness = 0;
		a.ta->next_generation[i].chromosomes = (int*) calloc(N, sizeof(int));
		a.ta->next_generation[i].index = i;
		a.ta->next_generation[i].chromosome_length = N;
		a.ta->next_generation[i].nr_chromosomes_1 = calloc(1, sizeof(int));
		a.ta->next_generation[i].nr_chromosomes_1[0] = 0;
	}

	pthread_barrier_wait(a.ta->barrier);////////////////////////////////

	// iterate for each generation
	for (int k = 0; k < a.ta->generations_count; ++k) {
		
		cursor = 0;
		
		compute_fitness_function_threaded(a.ta->objects, a.ta->current_generation, N, a.ta->sack_capacity, ID, P);

		pthread_barrier_wait(a.ta->barrier);////////////////////////////////

		// compute fitness and sort by it
		if (ID == 0) {
			qsort(a.ta->current_generation, N, sizeof(individual), cmpfunc);
		}

		pthread_barrier_wait(a.ta->barrier);//////////////////////////////////

		// keep first 30% children (elite children selection)
		count = N * 3 / 10;
		start = ID * (double)count / P;
		end = min((ID + 1) * (double)count / P, count);

		for (int i = start; i < end; ++i) {
			copy_individual(a.ta->current_generation + i, a.ta->next_generation + i);
		}
		cursor = count;

		// mutate first 20% children with the first version of bit string mutation
		count = N * 2 / 10;
		start = ID * (double)count / P;
		end = min((ID + 1) * (double)count / P, count);

		for (int i = start; i < end; ++i) {
			copy_individual(a.ta->current_generation + i, a.ta->next_generation + cursor + i);
			mutate_bit_string_1(a.ta->next_generation + cursor + i, k);
		}
		cursor += count;

		// mutate next 20% children with the second version of bit string mutation
		for (int i = start; i < end; ++i) {
			copy_individual(a.ta->current_generation + i + count, a.ta->next_generation + cursor + i);
			mutate_bit_string_2(a.ta->next_generation + cursor + i, k);
		}
		cursor += count;

		// crossover first 30% parents with one-point crossover
		// (if there is an odd number of parents, the last one is kept as such)
		count = N * 3 / 10;

		pthread_barrier_wait(a.ta->barrier);////////////////////////////////

		if (ID == 0) {
			if (count % 2 == 1) {
				copy_individual(a.ta->current_generation + N - 1, a.ta->next_generation + cursor + count - 1);
			}
		}

		pthread_barrier_wait(a.ta->barrier);////////////////////////////////

		count /= 2;
		start = ID * (double)(count / P);
		end = min((ID + 1) * (double)count / P, count);

		for (int i = start * 2; i < end * 2; i += 2) {
			crossover(a.ta->current_generation + i, a.ta->next_generation + cursor + i, k);
		}

		pthread_barrier_wait(a.ta->barrier);//////////////////////////////////

		// switch to new generation
		if (ID == 0) {
			tmp = a.ta->current_generation;
			a.ta->current_generation = a.ta->next_generation;
			a.ta->next_generation = tmp;
		}

		start = ID * (double)N / P;
		end = min((ID + 1) * (double)N / P, N);

		pthread_barrier_wait(a.ta->barrier);////////////////////////////////

		for (int i = start; i < end; ++i) {
			a.ta->current_generation[i].index = i;
		}

		if (ID == 0 && k % 5 == 0) {
			print_best_fitness(a.ta->current_generation);
		}

		pthread_barrier_wait(a.ta->barrier);////////////////////////////////

	}

	compute_fitness_function_threaded(a.ta->objects, a.ta->current_generation, N, a.ta->sack_capacity, ID, P);
	if (ID == 0) {
		qsort(a.ta->current_generation, N, sizeof(individual), cmpfunc);
		print_best_fitness(a.ta->current_generation);
	}

	pthread_barrier_wait(a.ta->barrier);////////////////////////////////

	free_generation_threaded(a.ta->current_generation, ID, P);
	free_generation_threaded(a.ta->next_generation, ID, P);

	return(0);
}

void run_genetic_algorithm(const sack_object *objects, int object_count, int generations_count, int sack_capacity, int nr_threads)
{
	pthread_barrier_t barrier;
	void* status;
	individual *current_generation = (individual*) calloc(object_count, sizeof(individual));
	individual *next_generation = (individual*) calloc(object_count, sizeof(individual));

	
	thread_argument ta;
	ta.barrier = &barrier;
	ta.current_generation = current_generation;
	ta.next_generation = next_generation;
	ta.objects = objects;
	ta.object_count = object_count;
	ta.generations_count = generations_count;
	ta.sack_capacity = sack_capacity;
	ta.nr_threads = nr_threads;

	thread_arg_with_id *a = (thread_arg_with_id *) calloc(nr_threads, sizeof(thread_arg_with_id));

	pthread_barrier_init(&barrier, NULL, nr_threads);

	int r;
	pthread_t threads[nr_threads];
	for (int i = 0; i < nr_threads; i++) {
		a[i].thread_id = i;
		a[i].ta = &ta;

        r = pthread_create(&threads[i], NULL, f, (void *)&a[i]);
		if (r) {
            printf("Eroare la crearea thread-ului %d\n", i);
            exit(-1);
		}
    }

	for (long id = 0; id < nr_threads; id++) {
			r = pthread_join(threads[id], &status);

			if (r) {
					printf("Eroare la asteptarea thread-ului %ld\n", id);
					exit(-1);
			}
  	}

	// free resources
	free(current_generation);
	free(next_generation);

	free(a);

	pthread_barrier_destroy(&barrier);
}