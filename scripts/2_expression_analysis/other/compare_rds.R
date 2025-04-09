data_list <- list.files(path = here("results/r_objects"),
                        pattern = ".RDS",
                        recursive = TRUE, full.names = TRUE) %>%
  purrr::set_names(., stringr::str_remove(basename(.), ".RDS")) %>%
  purrr::imap(~ readRDS(file = .x ))



data_list2 <- list.files(path = here("results/r_objects2"),
                        pattern = ".RDS",
                        recursive = TRUE, full.names = TRUE) %>%
  purrr::set_names(., stringr::str_remove(basename(.), ".RDS")) %>%
  purrr::imap(~ readRDS(file = .x ))

library(waldo)
waldo::compare(data_list, data_list2)
