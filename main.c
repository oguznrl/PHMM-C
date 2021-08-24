#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

enum TYPE {Beginning, Match, Insertion, Deletion, Unstable, Ending};

int ROW = 0, COLUMN = 0, UNSTABLE_COUNT = 0;
char *FILENAME = "data.txt";
/*---------------- PHMM Modelling Section  ----------------*/
struct state
{
    enum TYPE type;
    float probabilities[3];
    float row_lenght;
};

struct dna_sequence {
    char *sequence[30];
    char *path[30];
};

bool isNucleo(char *nucleo){
    return ( strcmp(nucleo, "A") == 0 ||
             strcmp(nucleo, "T") == 0 ||
             strcmp(nucleo, "G") == 0 ||
             strcmp(nucleo, "C") == 0
    );
}

int include(int *unstable_indexes,int column){
    for (int i = 0; i <UNSTABLE_COUNT ; ++i) {

        if(unstable_indexes[i]==column){

            return 1;
        }
    }
    return 0;
}

void get_shape(FILE *fp){
    char ch;
    while(!feof(fp)){
        ch=fgetc(fp);
        if(ch!=',') {
            COLUMN++;
        }
        if(ch=='\n'){
            ROW++;
        }
    }
    ROW++;
    COLUMN = (COLUMN-ROW) / ROW;
}

int * get_unstable_indexes(char *filename){
    FILE *fp;
    //reset
    UNSTABLE_COUNT = 0;
    char cr;
    int i,j;
    int *unstable_count = (int*) malloc(COLUMN * sizeof(int));
    int *unstable_row_indexes = (int*) malloc(COLUMN * sizeof(int));
    for (i = 0; i < COLUMN; i++) {
        unstable_count[i] = 0;
    }
    fp = fopen(filename, "r");
    while (!feof(fp)) {
        for (i = 0; i < ROW; i++) {
            for (j = 0; j < COLUMN; j++) {
                cr = getc(fp);
                if (cr == ',' || cr == '\n') {
                    j--;
                }
                else if (cr == '_') {
                    unstable_count[j] += 1;
                }
            }
        }
    }
    j = 0;
    for (i = 0; i < COLUMN; i++) {
        if(unstable_count[i] > 2){
            unstable_row_indexes[j] = i;
            UNSTABLE_COUNT+=1;
            j++;
        }
    }
    fclose(fp);

    return unstable_row_indexes;
}

void get_data(char *filename, struct dna_sequence *sequence) {
    FILE *fp;
    fp = fopen(filename, "r");
    char line[150];
    char *nucleotide;
    int _row = 0, _col = 0;

    for (_row = 0; _row < ROW; ++_row) {
        _col = 0;
        fgets(line, sizeof(line), fp);
        //printf("\nline : %s", line);

        nucleotide = strtok(line, ",");
        if(strcmp(nucleotide, "\n") != 0){
            sequence[_row].sequence[_col] =  strdup(nucleotide);
            //data[_row][_col] = strdup(nucleotide);
            //printf("nucleotid : %s | data : %s | col :%d\n", nucleotid, data[_row][_col], _col);
            _col++;
        }
        nucleotide = strtok(NULL, ",");

        while (nucleotide != NULL &&  strcmp(nucleotide, "\n") != 0) {
            sequence[_row].sequence[_col] =  strdup(nucleotide);
            //data[_row][_col] = strdup(nucleotide);
            //printf("nucleotid : %s | data : %s | col :%d\n", nucleotid, data[_row][_col], _col);
            _col++;
            nucleotide = strtok(NULL, ",");
        }

    }
    fclose(fp);
}

int len(struct state *state){
    int x=0;
    for (int i = 0; i <COLUMN ; ++i) {
        if(state[i].type==1||state[i].type==2||state[i].type==3)
            x++;
    }
    return x;
}

void create_model(struct state *match_states, struct state *insertion_states, struct state *deletion_states,struct state *beginnig_ending, struct dna_sequence *sequence){
    int *unstable_indexes = get_unstable_indexes(FILENAME);
    int m_d_count=0;
    int insertion_count=0;

    //Beginning inserstion 0
    struct state beginning=*(struct state*)malloc(sizeof(struct state));
    beginning.type=Beginning;
    beginning.probabilities[0]=0;beginning.probabilities[1]=0;beginning.probabilities[2]=0;
    beginning.row_lenght=0;
    beginnig_ending[0]=beginning;
    struct state insertion_state = *(struct state*) malloc(sizeof(struct state));
    insertion_state.type=Insertion;
    insertion_state.probabilities[0]=0;insertion_state.probabilities[1]=0;insertion_state.probabilities[2]=0;
    insertion_state.row_lenght=0;
    insertion_states[insertion_count]=insertion_state;
    insertion_count++;

    //ending

    struct state ending=*(struct state*)malloc(sizeof(struct state));
    ending.type=Ending;
    ending.probabilities[0]=0;ending.probabilities[1]=0;ending.probabilities[2]=0;
    ending.row_lenght=0;
    beginnig_ending[1]=ending;

    for (int index = 0; index <COLUMN; ++index) {
        if(!include(unstable_indexes,index)){
            int *row1=malloc(sizeof(int)*ROW);
            struct state match_state = *(struct state*) malloc(sizeof(struct state));
            match_state.type=Match;
            match_state.probabilities[0]=0;match_state.probabilities[1]=0;match_state.probabilities[2]=0;
            match_state.row_lenght=0;
            match_states[m_d_count]=match_state;
            struct state insertion_state = *(struct state*) malloc(sizeof(struct state));
            insertion_state.type=Insertion;
            insertion_state.probabilities[0]=0;insertion_state.probabilities[1]=0;insertion_state.probabilities[2]=0;
            insertion_state.row_lenght=0;
            insertion_states[insertion_count]=insertion_state;
            int *row3=malloc(sizeof(int)*ROW);
            struct state deletion_state = *(struct state*) malloc(sizeof(struct state));
            deletion_state.type=Deletion;
            deletion_state.probabilities[0]=0;deletion_state.probabilities[1]=0;deletion_state.probabilities[2]=0;
            deletion_state.row_lenght=0;
            deletion_states[m_d_count]=deletion_state;
            m_d_count++;
            insertion_count++;
        }

    }

    for (int _row = 0; _row < ROW ; ++_row) {
        struct dna_sequence _sequence = *(struct dna_sequence*) malloc(sizeof(struct dna_sequence));
        sequence[_row] = _sequence;
    }
}


/**
   * 1 - Elemanı al ve bir sonraki ile karşılaştır
   * 2 - Elemanı aldığın indeksin bulunduğu kolon unstable kolon mu kontrol et
   * 3 - Kolon Unstable değilse ve element _ değilse match
   * 4 - Kolon Unstable değilse ve element _ ise deletion
   * 5 - Kolon Unstable ise ve element _ değilse insertion
   * 6 - Kolon Unstable ise ve element _ ise ignore
   * 7 - Ek olarak her bir gittiği rotayı tut
   *
   * NOT:  Insertionlar arası geçiş yok eğer unstable ise ve önceki eleman
   * insertiona gittiyse kendisine dönecektir.
   */
char ***calculate_path(struct dna_sequence *sequence){
    int *unstable_indexes = get_unstable_indexes(FILENAME);
    int step;
    char ***path=(char ***)(malloc(sizeof(char **)*ROW));
    //char path[ROW][COLUMN][20];
    for (int i = 0; i <ROW ; ++i) {
        path[i]=(char **)(malloc(sizeof(char **)*COLUMN));
        step=0;
        for (int j = 0; j <COLUMN; ++j) {

            if (include(unstable_indexes,j)){

                if(isNucleo(sequence[i].sequence[j])){
                    path[i][j]=(char *)(malloc(sizeof(char *)*20));
                    sprintf(path[i][j],"i,%d",step);
                    printf("%s ",path[i][j]);
                }
                else{
                    path[i][j]=(char *)(malloc(sizeof(char *)*20));
                    sprintf(path[i][j],"X");
                    printf("%s ",path[i][j]);
                }
            }
            else{
                if(isNucleo(sequence[i].sequence[j])){
                    path[i][j]=(char *)(malloc(sizeof(char *)*20));
                    sprintf(path[i][j],"m,%d",step);
                    printf("%s ",path[i][j]);
                    step++;
                    }
                else{
                    path[i][j]=(char *)(malloc(sizeof(char *)*20));
                    sprintf(path[i][j],"d,%d",step);
                    printf("%s ",path[i][j]);
                    step++;
                }
            }
        }
        printf("\n");
    }
    return path;
}

void calculate_probability(char*** path,struct state *match_states, struct state *insertion_states, struct state *deletion_states,struct state *beginning_ending){
    char temp[20]="b";
    for (int i = 0; i <ROW ; ++i) {
        strcpy(temp,"b");
        for (int j = 0; j <COLUMN ; ++j) {
          if(temp[0]=='b' && path[i][j][0]=='m'){
              beginning_ending[0].probabilities[0]+=1;
              beginning_ending[0].row_lenght+=1;
              strcpy(temp,path[i][j]);}
          else if (temp[0]=='b' && path[i][j][0]=='i'){
              beginning_ending[0].probabilities[1]+=1;
              beginning_ending[0].row_lenght+=1;
              strcpy(temp,path[i][j]);}
          else if (temp[0]=='b' && path[i][j][0]=='d'){
              beginning_ending[0].probabilities[2]+=1;
              beginning_ending[0].row_lenght+=1;
              strcpy(temp,path[i][j]);}
          else if (temp[0]=='m' && path[i][j][0]=='m'){
               char *s;
               s=strtok(temp,",");
               s=strtok(NULL,",");
               match_states[atoi(s)].probabilities[0]+=1;
               match_states[atoi(s)].row_lenght+=1;
               strcpy(temp,path[i][j]);}
          else if (temp[0]=='m' && path[i][j][0]=='i'){
              char *s;
              s=strtok(temp,",");
              s=strtok(NULL,",");
              match_states[atoi(s)].probabilities[1]+=1;
              match_states[atoi(s)].row_lenght+=1;
              strcpy(temp,path[i][j]);}
          else if (temp[0]=='m' && path[i][j][0]=='d'){
              char *s;
              s=strtok(temp,",");
              s=strtok(NULL,",");
              match_states[atoi(s)].probabilities[2]+=1;
              match_states[atoi(s)].row_lenght+=1;
              strcpy(temp,path[i][j]);}
          else if (temp[0]=='i' && path[i][j][0]=='m'){
              char *s;
              s=strtok(temp,",");
              s=strtok(NULL,",");
              insertion_states[atoi(s)].probabilities[0]+=1;
              insertion_states[atoi(s)].row_lenght+=1;
              strcpy(temp,path[i][j]);}
          else if (temp[0]=='i' && path[i][j][0]=='i'){
              char *s;
              s=strtok(temp,",");
              s=strtok(NULL,",");
              insertion_states[atoi(s)].probabilities[1]+=1;
              insertion_states[atoi(s)].row_lenght+=1;
              strcpy(temp,path[i][j]);}
          else if (temp[0]=='i' && path[i][j][0]=='d'){
              char *s;
              s=strtok(temp,",");
              s=strtok(NULL,",");
              insertion_states[atoi(s)].probabilities[2]+=1;
              insertion_states[atoi(s)].row_lenght+=1;
              strcpy(temp,path[i][j]);}
          else if (temp[0]=='d' && path[i][j][0]=='m'){
              char *s;
              s=strtok(temp,",");
              s=strtok(NULL,",");
              deletion_states[atoi(s)].probabilities[0]+=1;
              deletion_states[atoi(s)].row_lenght+=1;
              strcpy(temp,path[i][j]);}
          else if (temp[0]=='d' && path[i][j][0]=='i'){
              char *s;
              s=strtok(temp,",");
              s=strtok(NULL,",");
              deletion_states[atoi(s)].probabilities[1]+=1;
              deletion_states[atoi(s)].row_lenght+=1;
              strcpy(temp,path[i][j]);}
          else if (temp[0]=='d' && path[i][j][0]=='d'){
              char *s;
              s=strtok(temp,",");
              s=strtok(NULL,",");
              deletion_states[atoi(s)].probabilities[2]+=1;
              deletion_states[atoi(s)].row_lenght+=1;
              strcpy(temp,path[i][j]);}}}
    for (int l = 0; l <ROW ; ++l) {
       char *s;
       s=strtok(path[l][COLUMN-1],",");
       if(s[0]=='m'){
           s=strtok(NULL,",");
           match_states[atoi(s)].probabilities[0]+=1;
           match_states[atoi(s)].row_lenght+=1;}
       else if(s[0]=='i'){
           s=strtok(NULL,",");
           insertion_states[atoi(s)].probabilities[0]+=1;
           insertion_states[atoi(s)].row_lenght+=1;}
       else if(s[0]=='d'){
           s=strtok(NULL,",");
           deletion_states[atoi(s)].probabilities[0]+=1;
           deletion_states[atoi(s)].row_lenght+=1;}}
    printf("\nBeginningn State\n\n");
    printf("beginning match->%.2f insertion->%.2f deletion->%.2f\n",beginning_ending[0].probabilities[0]/=beginning_ending[0].row_lenght ,beginning_ending[0].probabilities[1]/=beginning_ending[0].row_lenght ,beginning_ending[0].probabilities[2]/=beginning_ending[0].row_lenght);
    printf("\nMatch States\n\n");
    for (int k = 0; k <len(match_states) ; ++k) {
        if(match_states[k].row_lenght!=0){
            match_states[k].probabilities[0]/=match_states[k].row_lenght;
            match_states[k].probabilities[1]/=match_states[k].row_lenght;
            match_states[k].probabilities[2]/=match_states[k].row_lenght;}
        printf("match %d match->%.2f insertion->%.2f deletion->%.2f\n",k+1,match_states[k].probabilities[0],match_states[k].probabilities[1],match_states[k].probabilities[2]);}
    printf("\nDeletion States\n\n");
    for (int k = 0; k <len(deletion_states) ; ++k) {
        if(deletion_states[k].row_lenght!=0){
        deletion_states[k].probabilities[0]/=deletion_states[k].row_lenght;
        deletion_states[k].probabilities[1]/=deletion_states[k].row_lenght;
        deletion_states[k].probabilities[2]/=deletion_states[k].row_lenght;}
        printf("deletion %d match->%.2f insertion->%.2f deletion->%.2f\n",k+1,deletion_states[k].probabilities[0],deletion_states[k].probabilities[1],deletion_states[k].probabilities[2]);}
    printf("\nInsertion States\n\n");
    for (int k = 0; k <len(insertion_states) ; ++k) {
        if(insertion_states[k].row_lenght!=0){
            insertion_states[k].probabilities[0]/=insertion_states[k].row_lenght;
            insertion_states[k].probabilities[1]/=insertion_states[k].row_lenght;
            insertion_states[k].probabilities[2]/=insertion_states[k].row_lenght;}
        printf("insertion %d match->%.2f insertion->%.2f deletion->%.2f\n",k+1,insertion_states[k].probabilities[0],insertion_states[k].probabilities[1],insertion_states[k].probabilities[2]);}}

int main(void)
{
    FILE *fp;
    if ((fp = fopen(FILENAME, "r")) == NULL){
        printf("An error occured.");
        exit(1);}
    else {
        get_shape(fp);
        struct state *match_states = (struct state*) malloc(sizeof(struct state) * COLUMN);
        struct state *insertion_states = (struct state*) malloc(sizeof(struct state) * COLUMN);
        struct state *deletion_states = (struct state*) malloc(sizeof(struct state) * COLUMN);
        struct state *beginnig_ending=(struct state*) malloc(sizeof(struct state) * COLUMN);
        struct dna_sequence *dnaSequence = (struct dna_sequence*) malloc(sizeof(struct dna_sequence) * ROW);
        char ***path;
        create_model(match_states, insertion_states, deletion_states, beginnig_ending,dnaSequence);
        get_data(FILENAME, dnaSequence);
        path=calculate_path(dnaSequence);
        calculate_probability(path,match_states, insertion_states, deletion_states,beginnig_ending);}
    return 0;}

/**
 * 1 - Dosyayı oku (done)
 * 2 - Unstable durumlara karar ver (done)
 * 3 - Column sayısına göre dinamik dizileri oluştur (malloc) (done)
 * 4 - Row sayısı kadar dinamik sıralılar dizisi oluştur (done)
 * 5 - Tüm sıralıları teker teker kontrol et ve path dizisine stateleri ekle (done)
 * 6 - Sıra sıra ihtimalleri bulabilmek için kontrol et
 */