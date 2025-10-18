# packed_memory_array
A packed memory array for maintaining a sorted contiguous array in $$O(n)$$ space under inserts and deletions, while maintaining a constant number of gaps ($$O(1)$$ unused spaces) and requiring $$O(log(n)^2)$$ moves per `insert`/`delete` and allowing $$O(log(n))$$ `find` using a modified binary search.
