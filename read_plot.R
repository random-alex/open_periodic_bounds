#Compare open and periodic boundary conrition

# prereq ------------------------------------------------------------------

require(rhdf5)
require(multidplyr)
require(tidyverse)

dir <- 'data/a lot of temp and spin cor/'



# functions ---------------------------------------------------------------

my_read_matrix_param <- function(fol) {
  data1 <- h5read(fol,'simulation')
  df <- tibble(pairs = data1$results$`Spin Correlations`$labels,
               mean = data1$results$`Spin Correlations`$mean$value,
               error = data1$results$`Spin Correlations`$mean$error,
               error_convergence = data1$results$`Spin Correlations`$mean$error_convergence) %>% 
    rowwise() %>% 
    mutate(pairs = str_replace_all(pairs,'\\( | \\)',''),
           pairs = str_replace_all(pairs,' \\-\\- ',',')) %>% 
    ungroup() %>% 
    separate( 'pairs',c('x1','y1', 'x2','y2')) %>% 
    filter(x1 == '0' & y1 == '0' & x2 == '0' & y2 == '1')
  return(df)
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
my_read <- function(dir,cl = 3,pattern = '.out.h5', flag = 'standart'){
  
  cluster <- create_cluster(cl)
  if (flag == 'standart'){
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
    
  }
  else{
    df_data <- tibble(dir = list.files(dir,pattern = pattern,full.names = T),
                      id = rep(1:cl,length.out = NROW(dir))) %>% 
      partition(id,cluster = cluster) %>% 
      cluster_library(c('rhdf5','tidyverse')) %>% 
      cluster_assign_value('my_read_matrix_param',my_read_matrix_param) %>% 
      mutate(results = map(dir,my_read_matrix_param))  %>% 
      collect() %>% 
      ungroup() %>% 
      select(-id) %>% 
      unnest(results)  
  }
  
  
  parallel::stopCluster(cluster)
  return(df_data)
}


# read files --------------------------------------------------------------


df_init_par <- tibble(dir = list.files(dir,pattern = '.out.h5',full.names = T)) %>% 
  mutate(init_par = map(dir,my_read_parameters)) %>% 
  unnest(init_par) %>% 
  filter(parameter %in% c('L','T','LATTICE')) %>% 
  spread(parameter,value)


df_data <- my_read(dir,flag = 'matrix')


# work with staggered magnetization ---------------------------------------


df_data1 <- df_data %>% 
  filter(type %in% c("Staggered Magnetization Density^2")) %>% 
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
  # mutate(value = value/(as.numeric(as.character(L))^2)) %>%
  ggplot(aes(`T`,value,col = LATTICE)) +
  geom_line() +
  geom_pointrange(aes(ymax = value + error,ymin = value - error)) +
  geom_point(aes(shape = error_convergence),size = 3) +
  ylab("Staggered Magnetization Density^2") +
  facet_grid(type ~ L,scales = 'free',labeller = label_both) +
  theme_bw()


ss <- df %>% 
  # filter(type == 'Staggered Magnetization^2') %>% 
  group_by(L,`T`,LATTICE) %>% 
  summarise(n = n())



# work with Spin correlations --------------------------------------------------------

df <- left_join(df_init_par,df_data) %>% 
  select(-dir) %>% 
  arrange(as.numeric(L),as.numeric(T)) %>% 
  mutate(L = as_factor(L),
         error_convergence = as.factor(error_convergence),
         type = 'Spin Correlations',
         value = -mean)

df %>%
  mutate(`T` = as.numeric(`T`)) %>% 
  ggplot(aes(`T`,value,col = LATTICE)) +
  geom_line() +
  geom_pointrange(aes(ymax = value + error,ymin = value - error)) +
  geom_point(aes(shape = error_convergence),size = 3) +
  facet_grid(L ~ type,scales = 'free') +
  ylab('-Значения коррелятора на ближайших узлах') +
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
  mutate(value = value/as.numeric(as.character(L))^2) %>%
  mutate(value = (3/2 - (1 -1/3*mean(value[type == 'Magnetization^4'])/mean(value[type == 'Magnetization^2']^2)))) %>%
  ggplot(aes(`T`,value,col = L)) +
  geom_line() +
  # geom_pointrange(aes(ymax = value + error,ymin = value - error)) +
  geom_point(aes(shape = error_convergence),size = 3) +
  # facet_grid(type ~.,scales = 'free') +
  theme_bw()


