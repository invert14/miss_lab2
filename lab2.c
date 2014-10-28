#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

#define N 16771116

int fraction(mpz_t out, mpz_t num, mpz_t den, int d)
{
	mpz_t r;
	int zeros = 0;
	char beg = 1;
	mpz_init(r);
	mpz_set(r, num);
	mpz_set_ui(out, 0);
	while (d--)
	{
		mpz_mul_ui(num, r, 10);
		mpz_tdiv_qr(num, r, num, den);
		if (beg)
		{
			if (mpz_cmp_ui(num, 0) == 0)
				zeros++;
			else
				beg = 0;
		}
		mpz_abs(r, r);
		mpz_add(out, out, num);
		if (d)
			mpz_mul_ui(out, out, 10);
	}
	while (mpz_mod_ui(r, out, 10) == 0)
		mpz_divexact_ui(out, out, 10);
	return zeros;
}

void divAndPrint(mpz_t num, mpz_t den, int d)
{
	int zeros;
	char minus = 0;
	mpz_t rest, num_cp;
	mpz_init(num_cp);
	mpz_init(rest);
	mpz_tdiv_qr(num_cp, rest, num, den);
	if(mpz_cmp_ui(num_cp, 0) == 0 && mpz_cmp_ui(rest, 0) < 0)
		printf("-");
	mpz_abs(rest, rest);
	mpz_out_str(stdout, 10, num_cp);
	if(mpz_cmp_ui(rest, 0) > 0)
	{
		zeros = fraction(num_cp, rest, den, d);
		if(mpz_cmp_ui(num_cp, 0) > 0)
		{
			printf(".");
			while(zeros-- > 0) 
				printf("0");
			mpz_out_str(stdout, 10, num_cp);
		}
	}
	printf("\n");
	mpz_clear(rest);
	mpz_clear(num_cp);
}

/*****   macros create functional code   *****/
#define pivot_index() (begin+(end-begin)/2)
#define swap(a,b,t) ((t)=(a),(a)=(b),(b)=(t))
 
void sort(mpq_t array[], int begin, int end) {
   /*** Use of static here will reduce memory footprint, but will make it thread-unsafe ***/
	mpq_t pivot;
	mpq_init(pivot);
   //static int pivot;
   //static int t;     /* temporary variable for swap */

   if (end > begin) {
      int l = begin + 1;
      int r = end;
	mpq_swap(array[begin], array[pivot_index()]);
      //swap(array[begin], array[pivot_index()], t); /*** choose arbitrary pivot ***/
      mpq_set(pivot, array[begin]);
      while(l < r) {
	if (mpq_cmp(array[l], pivot) <= 0){
         //if (array[l] <= pivot) {
            l++;
         } else {
		while(l < --r && mpq_cmp(array[r], pivot) >= 0)
            //while(l < --r && array[r] >= pivot) /*** skip superfluous swaps ***/
               ;
            //swap(array[l], array[r], t); 
		mpq_swap(array[l], array[r]);
         }
      }
      l--;
      //swap(array[begin], array[l], t);
	mpq_swap(array[begin], array[l]);
      sort(array, begin, l);
      sort(array, r, end);
   }
}
 
#undef swap
#undef pivot_index

int y[10];
mpq_t *pi;

int y2[10];
mpq_t pi2[10];
int n2 = 0;

void fill_pi()
{
	pi = (mpq_t*) malloc (10 * sizeof(mpq_t));
	int i=0;
	for (i=0; i!=10; i++)
	{
		int up = (i+1)*(i+1) - i*i;
		mpq_init(pi[i]);
		mpq_set_ui(pi[i], up, 100);
		mpq_canonicalize(pi[i]);
	}
}

void fill_pi2()
{	
	int i=0;
	int fac = 1, fac2 = 2;

	mpq_t f;
	mpq_init(f);

	mpq_t f2;
	mpq_init(f2);

	mpq_t sum;
	mpq_init(sum);

	mpq_t one;
	mpq_init(one);
	mpq_set_ui(one, 1, 1);

	for (i=0; i!=9; i++)
	{
		mpq_init(pi2[i]);
		mpq_set_ui(f, 1, fac);
		mpq_set_ui(f2, 1, fac2);
		mpq_sub(pi2[i], f, f2);

		mpq_add(sum, sum, pi2[i]);
		fac *= (i+2);
		fac2 *= (i+3);
	}
	mpq_init(pi2[9]);
	mpq_sub(pi2[9], one, sum);
}

void print_table(mpq_t tab[], int len)
{
	int in;
	for (in=0; in!=len; in++)
	{
		gmp_printf("%Qd\n", tab[in]);
	}
}

void print_decimal(mpq_t akt, int d)
{
	mpz_t num;
	mpz_init(num);

	mpz_t den;
	mpz_init(den);

	mpq_get_num(num, akt);
	mpq_get_den(den, akt);
	divAndPrint(num, den, d);
}

int get_interval(mpq_t l)
{
	mpq_t pom;
	mpq_init(pom);
	mpq_set_ui(pom, 1, 100);
	int i=1;

	while (mpq_cmp(l, pom) >= 0)
	{
		i++;
		mpq_set_ui(pom, i*i, 100);
		mpq_canonicalize(pom);
	}
	i--;
	return i;
}

void update_series(int series)
{
	n2 += series;
	series--;
	if (series > 9)
		series = 9;
	y2[series]++;
}

int main(int argc, char* argv[])
{
	fill_pi();
	fill_pi2();

	int d;
	unsigned long int i = 0, n, p = 1;
	if(argc == 1) return 0;
	d = atoi(argv[1]);

	mpq_t *in;
	in = (mpq_t*) malloc (N * sizeof(mpq_t));
	
	mpq_t kminus;
	mpq_init(kminus);
	
	mpq_t kplus;
	mpq_init(kplus);

	mpq_t f;
	mpq_init(f);

	mpq_t akt;
	mpq_init(akt);

	mpq_t pom;
	mpq_init(pom);

	mpq_t ten;
	mpq_init(ten);
	mpq_set_ui(ten, 10, 1);

	double pom2;
	int interval;

	int series = 0;
	char omit = 0;
	char first = 1;


	while(gmp_scanf("%Qd", in[i]) != EOF)
	{
		mpq_canonicalize(in[i]);

		if (first == 1 ||  mpq_cmp(in[i], in[i-1]) > 0)
		{
			series++;
			first = 0;
		}
		else
		{
			update_series(series);

			series = 0;
			first = 1;
		}
		interval = get_interval(in[i]);
		if (interval < 10)
			y[interval]++;
		i++;
	}
	n = i;
	update_series(series);

	mpq_t np;
	mpq_init(np);
	mpq_canonicalize(np);

	mpq_t up;
	mpq_init(up);

	mpq_t(chi);
	mpq_init(chi);

	for (i=0; i!=10; i++)
	{
		//printf("%d: %d    *  ", i, y[i]);
		//gmp_printf("%Qd\n", pi[i]);

		mpq_set_ui(up, y[i], 1);
		
		mpq_set_ui(np, n, 1);
		mpq_mul(np, np, pi[i]);

		mpq_sub(up, up, np);
		mpq_mul(up, up, up);
		mpq_div(up, up, np);

		//gmp_printf("%Qd\n", up);

		mpq_add(chi, chi, up);
	}

	print_decimal(chi, d);

	//printf("n2=%d\n", n2);
	
	mpq_t chi2;
	mpq_init(chi2);
	for (i=0; i!=10; i++)
	{
		//printf("%d: %d    *  ", i, y2[i]);
		//gmp_printf("%Qd\n", pi2[i]);

		mpq_set_ui(up, y2[i], 1);
		
		mpq_set_ui(np, n2, 1);
		mpq_mul(np, np, pi2[i]);

		mpq_sub(up, up, np);
		mpq_mul(up, up, up);
		mpq_div(up, up, np);

		//gmp_printf("%Qd\n", up);

		mpq_add(chi2, chi2, up);
	}

	sort(in, 0, n);
	//print_table(in, n);

	mpq_t n_m;
	mpq_init(n_m);
	mpq_set_ui(n_m, n, 1);

	for (i=0; i!=n; i++)
	{
		mpq_set_ui(f, i, 1);
		mpq_div(f, f, n_m);
		mpq_sub(akt, in[i], f);
		
		if (mpq_cmp(akt, kminus) > 0)
			mpq_set(kminus, akt);
		
		//gmp_printf("%Qd - %Qd = %Qd = ", in[i], f, akt);
		//print_decimal(akt, d);
		
		mpq_set_ui(f, i+1, 1);
		mpq_div(f, f, n_m);
		mpq_sub(akt, f, in[i]);
		
		if (mpq_cmp(akt, kplus) > 0)
			mpq_set(kplus, akt);

		//gmp_printf("%Qd - %Qd = %Qd = ", f, in[i], akt);
		//print_decimal(akt, d);
	}

	print_decimal(kplus, d);
	print_decimal(kminus, d);
	
	print_decimal(chi2, d);


	return 0;
}
