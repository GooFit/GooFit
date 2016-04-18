
/*
 * strided_iterator.h
 *
 * original code: https://github.com/thrust/thrust/blob/master/examples/strided_range.cu
 *
 * Created on : Feb 25, 2016
 *      Author: Antonio Augusto Alves Junior
 */

/*
 *   This file is part of MCBooster.
 *
 *   MCBooster is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   MCBooster is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with MCBooster.  If not, see <http://www.gnu.org/licenses/>.
 */


#ifndef STRIDED_ITERATOR_H_
#define STRIDED_ITERATOR_H_

#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/iterator/permutation_iterator.h>
#include <thrust/functional.h>

namespace MCBooster{
/** \class strided_range
 * Strided range iterator original code: https://github.com/thrust/thrust/blob/master/examples/strided_range.cu
 */
template<typename Iterator>
class strided_range
{
public:

	typedef typename thrust::iterator_difference<Iterator>::type difference_type;

	struct stride_functor: public thrust::unary_function<difference_type,
			difference_type>
	{
		difference_type stride;

		stride_functor(difference_type stride) :
				stride(stride)
		{
		}

		__host__      __device__ difference_type operator()(
				const difference_type& i) const
		{
			return stride * i;
		}
	};

	typedef typename thrust::counting_iterator<difference_type> CountingIterator;
	typedef typename thrust::transform_iterator<stride_functor, CountingIterator> TransformIterator;
	typedef typename thrust::permutation_iterator<Iterator, TransformIterator> PermutationIterator;

	/// type of the strided_range iterator
	typedef PermutationIterator iterator;

	/// construct strided_range for the range [first,last)
	strided_range(Iterator first, Iterator last, difference_type stride) :
			first(first), last(last), stride(stride)
	{
	}

	iterator begin(void) const
	{
		return PermutationIterator(first,
				TransformIterator(CountingIterator(0), stride_functor(stride)));
	}

	iterator end(void) const
	{
		return begin() + ((last - first) + (stride - 1)) / stride;
	}

protected:
	Iterator first;
	Iterator last;
	difference_type stride;
};
}
#endif /* STRIDED_ITERATOR_H_ */
