### panaroo core edited

- [ ] Editar el script para que cuando de el roary.tsv para que recalcule las estadisticas para esas muestras.

- [ ] comprobar que en el roary no quita los duplicados.

- [x] change all_sample name to include core value

- [ ] add option to create tree of both, default core and personalised core

- [ ] add some fasta clean step (Fixcore?)

- [x] ~~test hmmclener in the server using a small aln: \~/perl5/bin/HmmCleaner.pl aln –~~

  - [ ] hmmer is its working on the server! But not in my computer.

- [ ] En teoria deberia de poner eliminar panaroo_edit core porque ya he creado el prev que se basa en que haya un panaroo previo

- hrnia que cambair los modules de baklta y shovill (al menos) para que sean modules individuales qu eno depended de mocules anteriores y utilizan una lista de accessions. Tengo que elegir si lista con path o lista de ficheros y directorio de entrada- tengo que ocmprobar si hay otros modulos que cmabiar. Pir ahora lo que la ia ha hehco e sun aprche para que utilice tanto los ahrdcoded o la lista de genomas dependienod del config, pero quiero que solo sea lista de gneomas.

- [ ] cambiar script de dowload para que el downoload de sra no lo haga en paralelo, solo 1 core, pero el fastq-dump si.