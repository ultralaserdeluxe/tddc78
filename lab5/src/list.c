#include <stdio.h>
#include <stdlib.h>

#include "coordinate.h" 

struct node {
  pcord_t* particle;
  struct node* prev;
  struct node* next;
};

struct node* create_list(){
  struct node* list = malloc(sizeof(struct node));
  list->particle = NULL;
  list->prev = NULL;
  list->next = NULL;
  return list;
}

void destroy_list(struct node* list){
  if(list->next != NULL) destroy_list(list->next);
  list->next = NULL;
  free(list);
}

int append_particle(pcord_t* particle, struct node* list){
  int length = 0;
  struct node* current = list;

  while(current->particle != NULL){
    current = current->next;
    length++;
  }

  current->particle = particle;
  current->next = create_list();

  return ++length;
}

void remove_particle(pcord_t* particle, struct node* list){
  struct node* current = list;

  while(current->particle != particle){
    current = current->next;
  }

  current->prev->next = current->next;
  current->next->prev = current->prev;

  free(current);
}

void print_list(struct node* list){
  struct node* current = list;

  printf("START->");
  do{
    if(current->particle != NULL) printf("(%f %f)->", current->particle->x, current->particle->y);
    current = current->next;
  }while(current != NULL);
  printf("END\n");
}

struct node* next(struct node* list){
  return list->next;
}

/* int main(){ */
/*   pcord_t p1  = {0, 1, 2, 3}; */
/*   pcord_t p2  = {1, 2, 3, 4}; */
/*   pcord_t p3  = {2, 3, 4, 5}; */

/*   struct node* list = create_list(); */
/*   print_list(list); */

/*   int length = append_particle(&p1, list); */
/*   printf("%d ", length); */
/*   print_list(list); */
  
/*   length = append_particle(&p2, list); */
/*   printf("%d ", length); */
/*   print_list(list); */

/*   length = append_particle(&p3, list); */
/*   printf("%d ", length); */
/*   print_list(list); */

/*   destroy_list(list); */
/* } */

