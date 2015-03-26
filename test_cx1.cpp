#include "cx1.h"

struct G {
	int a;
};

int main() {
	CX1<G, 65536> cx1;
	cx1.run();
	return 0;
}