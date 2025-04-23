#include <iostream>   
#include <vector>   
#include <gmp.h>   
#include <omp.h> 
#include <chrono>  
#include <cstdio>
using namespace std;
using namespace chrono;   
   
mpz_t k, step, step2, mult, mult2;   
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
bool TestEasy(vector <int>& prost, mpz_t Summ, mpz_t Razn) {  
    mpz_t del; 
    mpz_init(del);
    for (int i=0; i<prost.size(); i++) { // ! 
        mpz_set_ui(del, prost[i]); 
        if (mpz_divisible_p(Summ, del)){
            return false;
        }  
        if (mpz_divisible_p(Razn, del)){
            return false;
        } 
    }   
    mpz_clear(del);
    return true;
} 
 
bool TestMill(mpz_t Summ, mpz_t Razn) {
//summ = mult * 2^step + 1 
//берем случайное число 'a' от 998 до summ - 2, x = a^mult mod summ, если x == 1 или x == n - 1 число возможно простое, изменить 'a' 
//иначе x = x^2 mod summ, если x == 1, то число составное, если x==n-1 - число простое, количество итераций должно быть step-1 
//даже если x==n-1 при итерациях, надо продолжать вычисления с тем же 'a' 
mpz_t i, temp, ssum, sraz; 
mpz_init(temp);//темп это переменная - промежуток 
mpz_init(ssum); 
mpz_init(sraz); 
mpz_sub_ui(ssum, Summ, 1); //smul = summ - 1 
mpz_sub_ui(sraz, Razn, 1); //sraz = razn - 1 
mpz_sub_ui(temp, Summ, 1002); //summ - 2 - 998 - 2, позволит сразу взять расстояние от 998 до summ-2 
gmp_randstate_t state; 
gmp_randinit_mt(state);
gmp_randseed_ui(state, time(NULL));  
mpz_init(i);
mpz_t a, x, x2;//реализация для 6k+1, потом для 6k-1 
mpz_init(a); 
mpz_init(x); // 9 - 1
mpz_init(x2);
bool flag1 = false, flag2 = false;
//mpz_init_set_ui(dva, 2);
for (int j = 0; j <= 20; j++){ //основания для проверки  
    mpz_urandomm(a, state, temp); // создание рандомного a в промежутке от 998 до summ-4 
    mpz_add_ui(a, a, 998); 
    #pragma omp parallel
    {
    #pragma omp sections
    {
        #pragma omp section
        {
            mpz_powm(x, a, mult, Summ); // x = a^mult mod Summ  !
        }
        
        #pragma omp section
        {
            mpz_powm(x2, a, mult2, Razn); // x2 = a^mult2 mod Razn  !
        }
    }
    }
    if (mpz_cmp_ui(x, 1) == 0 || mpz_cmp(x, ssum) == 0) flag1 = true;
    if (mpz_cmp_ui(x2, 1) == 0 || mpz_cmp(x2, sraz) == 0) flag2 = true; 
    if (flag1 && flag2) return true;
    for (mpz_set_ui(i,1); mpz_cmp(i, step) < 0; mpz_add_ui(i, i, 1)) {//цикл для итераций, работает до step-1 раз 
    #pragma omp parallel
    {
    #pragma omp sections
    {
        #pragma omp section
        {
            mpz_mul(x, x, x);  // x=x^2
            mpz_fdiv_r(x, x, Summ); // x=x^2 mod summ   !
        }
        
        #pragma omp section
        {
            mpz_mul(x2, x2, x2);  
            mpz_fdiv_r(x2, x2, Razn); // x2=x^2 mod razn   !
        }
    }
}
        if(mpz_cmp(x, ssum) == 0) flag1 = true;
        if(mpz_cmp(x2, sraz) == 0) flag2 = true;
        if (flag1 && flag2) return true;
        if ((mpz_cmp_ui(x, 1) == 0)||(mpz_cmp_ui(x2, 1) == 0)) { 
            //gmp_fprintf(stdout, "Return false if a = %Zd, x = %Zd, x2 = %Zd, iterazia - %Zd \n", a, x, x2, i);
            return false; //число составное 
        } 
}    
}  
mpz_clear(a); 
mpz_clear(x); 
mpz_clear(x2);  
mpz_clear(i);
gmp_randclear(state);
mpz_clear(ssum); 
mpz_clear(sraz); 
mpz_clear(temp);
return true;
}
 
void RazlozSumm(){ // summ = mult * 2^step + 1   (16)17 - (18)19
    bool break1 = true;
    while (break1){ //mult%2==0 
        mpz_add_ui(step, step, 1); 
        mpz_divexact_ui(mult, mult, 2); 
        if (mpz_even_p(mult) == 0) {
            break1 = false;
        }
    } 
}  
void RazlozRazn(){ // razn = mult2*2^step2 + 1 
    bool break2 = true;
    while (break2){ 
        mpz_add_ui(step2, step2, 1); 
        mpz_divexact_ui(mult2, mult2, 2); 
        if (mpz_even_p(mult2) == 0) {
            break2 = false;
        }
    } 
} 
   
int main() {
    /*<bool> isPrime(10001, true);  // Решето Эратосфена для чисел до 10000   
    vector<int> prost;   
    isPrime[0] = isPrime[1] = false;
    for (int i = 2; i * i <= 10000; i++) {   
        if (isPrime[i]) {   
            for (int j = i * i; j <= 10000; j +=i) {
                isPrime[j] = false;   
            }   
        }   
    }   
    for (int i = 2; i <= 10000; i++) {
        if (isPrime[i]) prost.push_back(i);
    }   
    isPrime.clear();*/

   //2996863034895 * 2 ^ 1290000 / 6
    mpz_init(k);
    //Start();
    mpz_set_str(k, "1000000", 10);
    mpz_t Summ, Razn, k2;
    mpz_init_set_ui(k2, 0);
    mpz_init(mult);
    mpz_init(mult2);
    mpz_init(Summ); //k*6+1
    mpz_init(Razn); //k*6-1
    FILE* file = fopen("out.txt", "w");
    if (!file){
        perror("File can't open");
    }
    while (true) {
        mpz_add_ui(k, k, 1); // k++  
        //mpz_add_ui(k2, k2, 1);
        mpz_mul_ui(mult, k, 6); //k*6 18
        mpz_add_ui(Summ, mult, 1); // Summ = k*6+1 19
        mpz_sub_ui(Razn, mult, 1); // Razn = k*6-1 17
        mpz_sub_ui(mult2, Razn, 1);// mult2 = k*6-2 16
        int res1 = mpz_probab_prime_p(Razn, 5);
        if(res1 >= 1){
        //if (TestEasy(prost, Summ, Razn)) {
            int res2 = mpz_probab_prime_p(Summ, 5);
            
            if (res2 >= 1) {
                mpz_init_set_ui(step, 0);
                RazlozSumm();
                //mpz_add_ui(k2, k2, 1);
                gmp_fprintf(file,"%Zd \n * 2 ^ %Zd +-1 - Prime\n", mult, step); //mult*2^step +-1
            }
            //if(mpz_cmp_ui(k2, 10) == 0) break; 
        } 
    }
    fclose(file);
}