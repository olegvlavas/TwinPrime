#include <vector>   
#include <gmp.h>   
#include <omp.h> 
#include <ctime>  
#include <cstdio> 
using namespace std;   
  
mpz_t k;
void Start(){
    mpz_t base, exp;
    mpz_init(base);
    mpz_set_str(base, "2996863034895", 10);
    mpz_init_set_ui(exp, 2);
    mpz_pow_ui(k, exp, 1290000); // 2 ^ 1290000
    mpz_mul(k, k, base); // 2996863034895 * 2 ^ 1290000
    mpz_set_ui(exp, 6); // 2996863034895 * 2 ^ 1290000 / 6
    mpz_div(k, k, exp);
    mpz_clear(base);
    mpz_clear(exp);
}
 
void RazlozSumm(mpz_t mult1, mpz_t step1, FILE *file){ // summ = mult * 2^step + 1 
    while (true){ //mult%2==0 
        if (mpz_even_p(mult1) == 0) {
            break;
        }
        mpz_add_ui(step1, step1, 1); 
        mpz_divexact_ui(mult1, mult1, 2); 
    } 
    gmp_fprintf(file,"%Zd \n * 2 ^ %Zd +-1 \n", mult1, step1); //mult*2^step +-1 
    mpz_clear(mult1);
    mpz_clear(step1);
}  
int main() {
   // 2996863034895 * 2 ^ 1290000 / 6
    mpz_init(k);
  
    Start();
    //mpz_set_ui(k, 1000000);
  
    //clock_t start = clock();
    FILE* file = fopen("out.txt", "w");
    FILE* file2 = fopen("score.txt", "w");
    if (!file || !file2){
        perror("File can't open");
        return 1;
    }
    int score = 0;
    while (true) { 
        #pragma omp parallel for   
        for(int i = 1; i <= 1008; i++){
        int res1, res2;
        mpz_t k2, mult1, Summ1, Razn1;
        mpz_init(mult1);
        mpz_init(Summ1); //k*6+1   
        mpz_init(Razn1); //k*6-1 
        mpz_init(k2);
        mpz_add_ui(k2, k, i); // k++   
        mpz_mul_ui(mult1, k2, 6); //k*6   
        mpz_add_ui(Summ1, mult1, 1); // Summ = k*6+1   
        mpz_sub_ui(Razn1, mult1, 1); // Razn = k*6-1 
                res1 = mpz_probab_prime_p(Razn1, 5);
        if(res1 >= 1){
                res2 = mpz_probab_prime_p(Summ1, 5);
        if (res2 >= 1){
            mpz_t step1;
            mpz_init_set_ui(step1, 0);
            RazlozSumm(mult1, step1, file);
            mpz_clear(step1);
        }
        }
        fprintf(file2,"%d\n", score + i);
        mpz_clear(mult1);
        mpz_clear(Summ1);    
        mpz_clear(Razn1); 
        mpz_clear(k2);
        }
        score+=1008;
        fprintf(file2,"+ тысяча\n"); //mult*2^step +-1 
        mpz_add_ui(k, k, 1008);
        }
        //clock_t end = clock();
        //double elapsed_time = double(end - start) / CLOCKS_PER_SEC;
        //fprintf(file, "Execution time: %d seconds\n", elapsed_time);
        fclose(file);
    } 
