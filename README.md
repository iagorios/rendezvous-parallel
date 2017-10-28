# rendezvous-parallel

Instruções para compilar o serial

```
gcc rendezvous.c -o rendezvous.o -lm -fopenmp
```

Instruções para compilar o paralelo

```
gcc rendezvous.c -o rendezvous.o -lm -fopenmp
```

Instruções para executar o programa

```
./rendezvous.o [numerodeposicaoesiniciais]
```


