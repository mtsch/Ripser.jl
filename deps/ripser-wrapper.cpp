#include "ripser/ripser.cpp"
#include "ripser-wrapper.hpp"
/*
 * Simple wrapper to Ripser.
 */

void printmat(index_t n, value_t **mat) {
    for (int i = 0; i < n; ++i) {
        std::cout << i << ": ";
        for (int j = 0; j < i; ++j) {
            std::cout << mat[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

compressed_lower_distance_matrix from_arrays(index_t n, value_t **mat) {
	std::vector<value_t> distances;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <= i; ++j) {
            distances.push_back(mat[i][j]);
        }
    }
	return compressed_lower_distance_matrix(std::move(distances));
}

// Essentially a copy of the main function that takes arguments.
void ripser(index_t len, value_t **mat, index_t dim_max, value_t threshold, coefficient_t modulus) {

    compressed_lower_distance_matrix dist = from_arrays(len, mat);

	index_t n = dist.size();

	std::cout << "distance matrix with " << n << " points" << std::endl;

	auto value_range = std::minmax_element(dist.distances.begin(), dist.distances.end());
	std::cout << "value range: [" << *value_range.first << "," << *value_range.second << "]" << std::endl;

	dim_max = std::min(dim_max, n - 2);

	binomial_coeff_table binomial_coeff(n, dim_max + 2);
	std::vector<coefficient_t> multiplicative_inverse(multiplicative_inverse_vector(modulus));

	std::vector<diameter_index_t> columns_to_reduce;

	{
		union_find dset(n);
		std::vector<diameter_index_t> edges;
		rips_filtration_comparator<decltype(dist)> comp(dist, 1, binomial_coeff);
		for (index_t index = binomial_coeff(n, 2); index-- > 0;) {
			value_t diameter = comp.diameter(index);
			if (diameter <= threshold) edges.push_back(std::make_pair(diameter, index));
		}
		std::sort(edges.rbegin(), edges.rend(), greater_diameter_or_smaller_index<diameter_index_t>());

#ifdef PRINT_PERSISTENCE_PAIRS
		std::cout << "persistence intervals in dim 0:" << std::endl;
#endif

		std::vector<index_t> vertices_of_edge(2);
		for (auto e : edges) {
			vertices_of_edge.clear();
			get_simplex_vertices(get_index(e), 1, n, binomial_coeff, std::back_inserter(vertices_of_edge));
			index_t u = dset.find(vertices_of_edge[0]), v = dset.find(vertices_of_edge[1]);

			if (u != v) {
#ifdef PRINT_PERSISTENCE_PAIRS
				if (get_diameter(e) > 0) std::cout << " [0," << get_diameter(e) << ")" << std::endl;
#endif
				dset.link(u, v);
			} else
				columns_to_reduce.push_back(e);
		}
		std::reverse(columns_to_reduce.begin(), columns_to_reduce.end());

#ifdef PRINT_PERSISTENCE_PAIRS
		for (index_t i = 0; i < n; ++i)
			if (dset.find(i) == i) std::cout << " [0, )" << std::endl << std::flush;
#endif
	}

	for (index_t dim = 1; dim <= dim_max; ++dim) {
		rips_filtration_comparator<decltype(dist)> comp(dist, dim + 1, binomial_coeff);
		rips_filtration_comparator<decltype(dist)> comp_prev(dist, dim, binomial_coeff);

		hash_map<index_t, index_t> pivot_column_index;
		pivot_column_index.reserve(columns_to_reduce.size());

		compute_pairs(columns_to_reduce, pivot_column_index, dim, n, threshold, modulus, multiplicative_inverse, dist,
		              comp, comp_prev, binomial_coeff);

		if (dim < dim_max) {
			assemble_columns_to_reduce(columns_to_reduce, pivot_column_index, comp, dim, n, threshold, binomial_coeff);
		}
	}
}
