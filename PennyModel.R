# Nbit - genom lenght
# N0   - start population size
# Nmax - max population size (Typical Nmax = 10 · N0)
# m    - the probability of additional mutations obtained at birth
# M    - The amount of additional mutations obtained at birth
# T    - the body's resistance to "diseases"
# Rmin - age of reaching sexual maturity
# Rmax - age at which an individual loses its ability to procreate (usually Rmax = Nbit)
# b    - the probability of reproduction
# B    - the number of offspring produced by an individual during the reproductive period
# p    - probability of bad->good mutation
# P    - number good mutations
# Tmax - simulation time

model <- function(Nbit, N0, Nmax, m, M, T, Rmin, Rmax, b, B, p, P, Tmax) {
  subjects <- matrix(as.integer(runif(N0 * (Nbit + 2)) < 0.1), nrow = N0, ncol = Nbit + 2)
  subjects[, Nbit + 1] <- 0
  subjects[, Nbit + 2] <- TRUE
  subjects <- data.frame(subjects)
  View(subjects)
  next_gen_id <- 1

  for (t in 1:Tmax) {
    subjects_next_gen <- subjects
    #if(nrow(subjects) > 50000) break

    subjects$X17 <- subjects$X17 + 1
    for (s in nrow(subjects)) {
      if (subjects[s, 'X17'] > Nbit) { #check that the maximum age is exceeded
        subjects[s, 'X18'] <- FALSE
      }
    }
    for (s in nrow(subjects)) {
      if (subjects[s, 'X18']) {
        if (sum(subjects[s, 1:subjects[s, 'X17']]) > T) subjects[s, 'X18'] <- FALSE #check that the number of exposed and active genes does not exceed body resistance
      }
    }
    subjects <- verhulst(subjects, Nmax) # Verhulst


    for (s in nrow(subjects)) {
      if (check_population_size(subjects, Nmax, Nbit)) {
        if ((subjects[s, 'X18']) &
          (subjects[s, 'X17'] >= Rmin) &
          (subjects[s, 'X17'] <= Rmax)) { #reproduction

          for (i in range(1, B)) {
            if (draw(b)) {
              subjects_next_gen[next_gen_id,] <- c(subjects[s, 1:Nbit], 0, 1)
            }else {
              #subjects_next_gen[next_gen_id,] <- c(subjects[s, 1:Nbit], 0, 1)
            }
            print("id : ", next_gen_id)
            print(subjects_next_gen[next_gen_id,])
            print(subjects[s,])
            next_gen_id <- next_gen_id + 1
          }
        }
      }
    }
    subjects <- to_the_next_gen(subjects, subjects_next_gen, Nbit, t)
    subjects <- mutacje(subjects, M, m, P, p, Nbit)
  }
  View(subjects)
  return(subjects)
}

is.not.null <- function(x) !is.null(x)

draw <- function(propability) {
  propability <- propability / 100
  if (propability > 1) propability <- 1
  smp <- sample(c(1, 0), size = 10, replace = TRUE, prob = c(propability, 1 - propability))
  if (sample(smp, 1) == 1) return(TRUE)
  else return(FALSE)
}

add_mutation <- function(parent_genome, m, M) {
  pos <- 1
  added_mutations <- 0
  output <- list()
  parent_genome <- as.numeric(parent_genome)
  for (g in 1:(length(parent_genome))) {
    if (added_mutations < M) {
      if (parent_genome[g] == 0) {
        if (draw(m)) {
          output[pos] <- 1
          added_mutations <- added_mutations + 1
        }else {
          output[pos] <- 0
        }
      }else {
        output[pos] <- 1
      }
    }
    pos <- pos + 1
  }
  return(output)
}

verhulst <- function(mtrx, Nmax) {
  for (s in nrow(mtrx)) {
    mtrx[s, 'X18'] <- !draw(nrow(mtrx) / Nmax)
  }
  return(mtrx)
}

check_population_size <- function(subjects, Nmax, Nbit) {
  dead <- 0
  for (s in nrow(subjects)) {
    if (!subjects[s, 'X18']) dead <- dead + 1
  }
  if (Nmax > nrow(subjects) - dead) return(TRUE)
  else return(FALSE)
}

mutacje <- function(mtrx, M, m, P, p, NBit) {
  for (i in nrow(mtrx)) {
    il <- sum(runif(M) < m)
    if (il > 0) mtrx[i, sample(NBit)[1:il]] <- 1
    il <- sum(runif(P) < p)
    if (il > 0)  mtrx[i, sample(NBit)[1:il]] <- 0
  }
  return(mtrx)
}

to_the_next_gen <- function(subject, next_gen, Nbit, t) {
  if (t == 1) {
    not_dead <- sum(subject$X18)
    if (not_dead > 0) {
      output_matrix <- subject[subject$X18 == TRUE,]
    }else {
      print('all subject died')
      sys.on.exit()
    }
  }else {
    not_dead <- sum(subject$X18) + sum(next_gen$X18)
    if (not_dead > 0) {
      out1 <- subject[subject$X18 == TRUE,]
      out2 <- next_gen[next_gen$X18 == TRUE,]

      output_matrix <- rbind(out1, out2)
    }else {
      print('all subject died')
      sys.on.exit()
    }
  }

  return(output_matrix)
}
make_graph <-function (model.res){
  len <- model.res$X17[1]
  plot_vec <- vector(length = len)
  for (i in 1:len){
    plot_vec[i] <- sum(model.res$X17 == i)
  }
  plot(plot_vec,xlab = 'time',ylab = 'population size',col = "red",type = 'b')
}

model.resolut <- model(16, 6, 60, 20, 10, 10, 10, 16, 50, 6, 20, 6, 15)

make_graph(model.resolut)
