# Simulazione_-_Controllo_Manipolatore_Planare_2_bracci
L'obiettivo è quello di sviluppare un applicativo in ambiente MATLAB che affronti la pianificazione, la simulazione e il controllo del robot.

Definizione del Modello Dinamico: Si definisce un sistema a due gradi di libertà, fornendo i parametri fisici necessari come le lunghezze dei link ($a_1 = 1 m, a_2 = 0.6 m$), le masse ($5 kg$ e $3 kg$), le inerzie e le posizioni dei centri di massa. Il progetto si basa sul modello dinamico completo, governato dall'equazione $B(q)\dot{q} + C(q,\dot{q})\dot{q} + g(q) = \tau$, le cui matrici per l'inerzia $B(q)$, per le forze di Coriolis $C(q,\dot{q})$ e per il vettore di gravità $g(q)$ sono già interamente modellate ed esplicitate nel testo (assumendo di trascurare inerzia dei motori e attriti).  

Fase 1 - Pianificazione delle Traiettorie: Il primo task pratico consiste nel generare delle traiettorie punto-punto per far muovere il braccio da una configurazione iniziale $q_0 = [\pi/4, \pi/3]$ verso una posizione finale in un lasso di tempo di 10 secondi. Le traiettorie devono seguire due leggi orarie distinte:  Per il giunto 1 bisogna raggiungere la coordinata $\theta_1 = \pi/2$ utilizzando un polinomio di interpolazione cubica.  Per il giunto 2 bisogna raggiungere la coordinata $\theta_2 = -\pi/2$ sfruttando invece un profilo di velocità trapezoidale.  

Fase 2 - Simulazione: Sulla base delle traiettorie pianificate e del modello matematico, è necessario implementare la simulazione del robot utilizzando la dinamica diretta.  

Fase 3 - Sistema di Controllo: Per permettere al robot di seguire le traiettorie calcolate e vincere le forze in gioco, viene richiesto di implementare un controllore PD (Proporzionale-Derivativo) arricchito da una compensazione di gravità nello spazio dei giunti.
