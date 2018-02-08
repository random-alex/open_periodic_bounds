#Compare open and periodic boundary conrition

# prereq ------------------------------------------------------------------

require(rhdf5)
require(multidplyr)
require(tidyverse)

dir <- 'data/'



# functions ---------------------------------------------------------------

my_prep_pairs <- function(pairs) {
  str_extract_all(pairs,'[0-9]',simplify = F) %>%
    unlist() %>% 
    .[c(3,4)] %>% 
    str_c(sep = '_',collapse = '_')
}
my_prep_plot <- function(df,t) {
  temp1 <- df %>% 
    filter(temp == t) %>% 
    .$mean
  dim(temp1) <- c(NROW(unique(df$x)),NROW(unique(df$y)))
  return(temp1)
}
my_read_parameters <- function(fol){
  h5read(fol,'parameters') %>% 
    enframe(name = 'parameter') %>%
    mutate(value = map(value,as.character)) %>% 
    unnest(value)
  
}
my_read_sim <- function(fol) {
  data1 <- h5read(fol,'simulation')$results %>% 
    enframe(name = 'type') %>% 
    mutate(value = map(value,enframe,name = 'parameter')) 
  return(data1)
}
my_read <- function(dir,cl = 3,pattern = '.out.h5'){
  
  cluster <- create_cluster(cl)
  
  df_data <- tibble(dir = list.files(dir,pattern = pattern,full.names = T),
                    id = rep(1:cl,length.out = NROW(dir))) %>% 
    partition(id,cluster = cluster) %>% 
    cluster_library(c('rhdf5','tidyverse')) %>% 
    cluster_assign_value('my_read_sim',my_read_sim) %>% 
    mutate(results = map(dir,my_read_sim))  %>% 
    collect() %>% 
    ungroup() %>% 
    select(-id) %>% 
    unnest(results)  
  parallel::stopCluster(cluster)
  return(df_data)
}



# read files --------------------------------------------------------------


df_init_par <- tibble(dir = list.files(dir,pattern = '.out.h5',full.names = T)) %>% 
  mutate(init_par = map(dir,my_read_parameters)) %>% 
  unnest(init_par) %>% 
  filter(parameter %in% c('L','T','LATTICE')) %>% 
  spread(parameter,value)


df_data <- my_read(dir)


# work with staggered magnetization ---------------------------------------


df_data1 <- df_data %>% 
  filter(type %in% c("Staggered Magnetization^2")) %>% 
  unnest(value) %>% 
  filter(parameter %in% c('mean')) %>% 
  select(-c(parameter) ) %>% 
  mutate(res = map(value,enframe)) %>% 
  unnest(res) %>% 
  unnest(value) %>% 
  spread(name,value)

df <- left_join(df_init_par,df_data1) %>% 
  select(-dir) %>% 
  arrange(as.numeric(L),as.numeric(`T`)) %>% 
  mutate(L = as_factor(L),
         error_convergence = as.factor(error_convergence))

df %>% 
  filter(error_convergence == 1)


df %>%
  mutate(`T` = as.numeric(`T`)) %>% 
  group_by(L,`T`) %>% 
  mutate(value = value/as.numeric(L)^2) %>%
  ggplot(aes(`T`,value,col = L)) +
  geom_line() +
  geom_pointrange(aes(ymax = value + error,ymin = value - error)) +
  geom_point(aes(shape = error_convergence),size = 3) +
  facet_grid(type ~ LATTICE,scales = 'free') +
  theme_bw()


ss <- df %>% 
  # filter(type == 'Staggered Magnetization^2') %>% 
  group_by(L,`T`,LATTICE) %>% 
  summarise(n = n())



# work with energy --------------------------------------------------------


df_data1 <- df_data %>% 
  filter(type %in% c('Energy')) %>% 
  unnest(value) %>% 
  filter(parameter %in% c('mean')) %>% 
  select(-c(parameter) ) %>% 
  mutate(res = map(value,enframe)) %>% 
  unnest(res) %>% 
  unnest(value) %>% 
  spread(name,value)

df <- left_join(df_init_par,df_data1) %>% 
  select(-dir) %>% 
  arrange(as.numeric(L),as.numeric(T)) %>% 
  mutate(L = as_factor(L),
         error_convergence = as.factor(error_convergence))

df %>%
  mutate(`T` = as.numeric(`T`)) %>% 
  group_by(L,`T`) %>% 
  ggplot(aes(`T`,value,col = LATTICE)) +
  geom_line() +
  geom_pointrange(aes(ymax = value + error,ymin = value - error)) +
  geom_point(aes(shape = error_convergence),size = 3) +
  facet_grid(L ~ type,scales = 'free') +
  theme_bw()



# work with Magnetization^2 -----------------------------------------------

df_data1 <- df_data %>% 
  filter(type %in% c('Magnetization^2','Magnetization^4')) %>% 
  unnest(value) %>% 
  filter(parameter %in% c('mean')) %>% 
  select(-c(parameter) ) %>% 
  mutate(res = map(value,enframe)) %>% 
  unnest(res) %>% 
  unnest(value) %>% 
  spread(name,value)

df <- left_join(df_init_par,df_data1) %>% 
  select(-dir) %>% 
  arrange(as.numeric(L),as.numeric(`T`)) %>% 
  mutate(L = as_factor(L),
         error_convergence = as.factor(error_convergence))


df %>% 
  filter(error_convergence == 1)


df %>%
  # filter(type == 'Staggered Magnetization^2') %>% 
  mutate(`T` = as.numeric(`T`)) %>% 
  group_by(L,`T`) %>% 
  mutate(value = (3/2 - (1 -1/3*mean(value[type == 'Magnetization^4'])/mean(value[type == 'Magnetization^2']^2)))) %>%
  ggplot(aes(`T`,value,col = L)) +
  geom_line() +
  # geom_pointrange(aes(ymax = value + error,ymin = value - error)) +
  geom_point(aes(shape = error_convergence),size = 3) +
  # facet_grid(type ~.,scales = 'free') +
  theme_bw()


