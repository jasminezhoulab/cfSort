#ifndef TEMPLATE_UTILS_H
#define TEMPLATE_UTILS_H

// print vector of any primary data type to the format, e.g., 10469,10471,10484,10489,10493,10497,10525,10542
template <typename T> void print_vec(ostream& of, vector<T> & v, string delimit=",", string prefix="") {
	switch (v.size()) {
		case 0:
			return;
		case 1:
			of << prefix << v[0];
			break;
		default: // v.size()>=2
			int i;
			of << prefix;
			for (i=0; i<v.size()-1; i++) of << v[i] << delimit;
			of << v[i];
	}
}

#endif

