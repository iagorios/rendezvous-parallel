# rendezvous-parallel

Preparando a execução

1- Clonar o repositório <br />
2- Inserir o arquivo de entrada com o nome "in.dat" na pasta clonada.

Executando o código serial
```
gcc rendezvous-serial.c -o rendezvous-serial.o -lm

./rendezvous-serial.o [numeroDePosicoesIniciais]

```

Executando o código paralelo
```
gcc rendezvous-openMP.c -o rendezvous-openMP.o -lm -fopenmp

./rendezvous-serial.o [numeroDePosicoesIniciais]
```

O parâmetro ```-pg``` deve ser utilizado para se obter a análise do gprof.
