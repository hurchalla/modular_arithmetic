MontySqrtRange is an optimization of montgomery arithmetic that applies only in the case where the modulus is less than the square root of R.  R represents the value 2^(bit_width_of_type_used) - for example if the type used is uint64_t, R = 2^64.

All files in this experimental directory are expected to be good and immediately usable if desired.  However, I consider them to be experiments that aren't useful enough to be moved outside of an "experimental" directory.
