# rendezvous-parallel

Instruções para compilar o serial

```
gcc rendezvous.c -o rendezvous.o -lm
```

Instruções para compilar o paralelo

```
gcc rendezvous-parallel.c -o rendezvous.o -lm -fopenmp
```

Instruções para executar o programa

```
./rendezvous.o [numerodeposicaoesiniciais]
```


