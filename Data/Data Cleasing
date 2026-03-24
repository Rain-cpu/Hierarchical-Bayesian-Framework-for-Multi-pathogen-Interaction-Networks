####Data Cleasing####
####By Suyi####
####20250411####


# Load required packages --------------------------------------------------
library(readxl)
library(dplyr)
library(stringr)

# 1. Read the Excel sheet ----------------------------------------------------
file_path <- "file_path"
df_raw <- read_excel(file_path, sheet = sheet)

# Step 1: Assign layer type
classify_layer <- function(index, unit) {
  if (index == "OR") {
    return("OR")
  } else if (!is.na(unit) && tolower(unit) == "log10") {
    return("LOG")
  } else {
    return("VIC")
  }
}

df <- df_raw %>%
  mutate(layer = mapply(classify_layer, `OR/Index`, unit))


#  2. Generate from-to pairs and determine directionality -----------------
parse_pair <- function(pair_str) {
  parts <- unlist(str_split(pair_str, "-"))
  if (length(parts) == 2) return(parts) else return(c(NA, NA))
}

pair_info <- t(sapply(df$`Virus-pair`, parse_pair))
colnames(pair_info) <- c("pair1", "pair2")
df$from <- ifelse(!is.na(df$From), df$From, pair_info[,1])
df$to   <- ifelse(!is.na(df$To), df$To, pair_info[,2])

# Add edge directionality indicator
# If both From and To are NA → undirected, otherwise directed
df$direction <- ifelse(is.na(df$From) & is.na(df$To), "undirected", "directed")



# 3. Compute or assign SE -------------------------------------------------

# SE from 95% CI (if lower/upper given)
compute_se_from_ci <- function(lower, upper) {
  return((log(upper) - log(lower)) / (2 * 1.96))
}

# Default SE if no CI is available【SE need sensitivity analysis】
assign_se_raw <- function(layer, se_input, lower, upper) {
  if (layer == "OR") {
    if (!is.na(lower) && !is.na(upper)) {
      return(compute_se_from_ci(lower, upper))
    } else {
      return(se_input)  # fallback to manual SE if CI not available
    }
  } else if (layer == "LOG") {
    return(1.5) #set log se as 1.5
  } else if (layer == "VIC") {
    return(1.0) #set vic se as 1.0
  } else {
    return(NA)
  }
}

# Apply SE assignment across rows
df$se <- NA
df$se <- mapply(assign_se_raw, df$layer, df$se, df$CI_lower, df$CI_upper)

df <- df %>%
   mutate(
    sigma = ifelse(!is.na(IF) & !is.na(se), se / sqrt(log(IF + 1)), NA)
    #sigma = ifelse(!is.na(IF) & !is.na(se), se, NA)
  ) %>%
  mutate(
    log_value = case_when(
      layer == "OR"  ~ log(Value),  # OR: log
      layer == "LOG" ~ Value,       # LOG: log10 as data unit
      layer == "VIC" ~ log(Value),  # VIC: log
      TRUE           ~ NA_real_
    )
  )

cleaned_data <- df %>%
  select(from, to, direction, layer, value = Value, log_value, se, sigma, IF = `IF`, source = Num)

#Standardization of log layer
log_range <- cleaned_data %>%
  filter(layer == "LOG") %>%
  summarise(
    min_neg = min(log_value[log_value < 0], na.rm = TRUE),  # e.g. -5.0
    max_pos = max(log_value[log_value > 0], na.rm = TRUE)   # e.g. +1.0
  )

#linear mapping
# [min_neg, 0] → [-2, 0]
# [0, max_pos] → [0, +0.3]
cleaned_data <- cleaned_data %>%
  mutate(
    log_value_scaled = case_when(
      layer != "LOG" ~ log_value,  
      log_value < 0  ~ (log_value - 0) / (log_range$min_neg - 0) * (-2),  
      log_value > 0  ~ (log_value / log_range$max_pos) * 0.3,
      log_value == 0 ~ 0,
      TRUE ~ NA_real_
    )
  )

# 5. Summary and inspection -----------------------------------------------
# Step 
#summary(cleaned_data)
#table(cleaned_data$layer)
#head(cleaned_data)

# Optional: Save cleaned data
# write.csv(cleaned_data, "prepared_vci_inputs.csv", row.names = FALSE)

# 6. Prepare Stan/INLA input format and flu_all constraint data -----------
# 6.1 Split data into network input and flu_all constraint input ----------
excluded_node <- "IV"
virus_nodes <- sort(unique(c(cleaned_data$from, cleaned_data$to)))
virus_nodes <- setdiff(virus_nodes, excluded_node)
pair_id <- data.frame(virus = virus_nodes, id = seq_along(virus_nodes))

# Split IV rows as separate flu_all observations
flu_all_data <- cleaned_data %>% filter(from == "IV" | to == "IV")
network_data <- cleaned_data %>% filter(!(from == "IV" | to == "IV"))

# Assign node IDs (for network_data only)
network_data <- network_data %>%
  left_join(pair_id, by = c("from" = "virus")) %>% rename(from_id = id) %>%
  left_join(pair_id, by = c("to" = "virus")) %>% rename(to_id = id)
cleaned_data <- cleaned_data %>%
  left_join(pair_id, by = c("from" = "virus")) %>% rename(from_id = id) %>%
  left_join(pair_id, by = c("to" = "virus")) %>% rename(to_id = id)


# 6.2 Assign pairwise ID for unique virus combinations --------------------
network_data <- network_data %>%
  mutate(
    pair_key = ifelse(direction == "undirected",
                      paste0(pmin(from_id, to_id), "_", pmax(from_id, to_id)),
                      paste0(from_id, "_", to_id))
  )

pair_lookup <- data.frame(pair_key = unique(network_data$pair_key)) %>%
  mutate(pair_index = row_number())

network_data <- network_data %>%
  left_join(pair_lookup, by = "pair_key")

# Remove any rows where from_id or to_id is NA (e.g. from IV being excluded)
network_data <- network_data %>% filter(!is.na(from_id) & !is.na(to_id))


# 6.3 Retain IV as covariate information ----------------------------------
network_data <- network_data %>% mutate(has_IV = 0)

# Final structure for Stan/INLA
stan_input <- network_data %>%
  mutate(
    layer_or = ifelse(layer == "OR", 1, 0),
    layer_log = ifelse(layer == "LOG", 1, 0),
    layer_vic = ifelse(layer == "VIC", 1, 0)
  ) %>%
  select(pair_index, from_id, to_id, direction, layer, layer_or, layer_log, layer_vic, has_IV, log_value_scaled, se, sigma, IF, source)


# 6.4 Prepare flu_all_data for constraint modeling (only OR layer with direction) --------
flu_all_or <- flu_all_data %>%
  filter(layer == "OR", direction == "directed", from == "IV" | to == "IV") %>%
  mutate(
    flu_partner = ifelse(from == "IV", to, from),                     
    direction_type = ifelse(from == "IV", "outgoing", "incoming"),   
    log_value_scaled = log(value),
    sigma = se / sqrt(log(IF + 1))                                    
    #sigma = se                                   #Sensitivity analysis
  ) %>%
  left_join(pair_id, by = c("flu_partner" = "virus")) %>%
  rename(flu_partner_id = id) %>%
  mutate(
    IAV_id = 2,
    IBV_id = 3,
    key_a = ifelse(direction_type == "outgoing", 
                   paste0(IAV_id, "_", flu_partner_id), 
                   paste0(flu_partner_id, "_", IAV_id)),
    key_b = ifelse(direction_type == "outgoing", 
                   paste0(IBV_id, "_", flu_partner_id), 
                   paste0(flu_partner_id, "_", IBV_id))
  ) %>%
  left_join(pair_lookup, by = c("key_a" = "pair_key")) %>%
  rename(flu_A_idx = pair_index) %>%
  left_join(pair_lookup, by = c("key_b" = "pair_key")) %>%
  rename(flu_B_idx = pair_index) %>%
  mutate(flu_is_directed = 1) %>%
  select(flu_partner, flu_A_idx, flu_B_idx, log_value_scaled, sigma, flu_is_directed)
#write.csv(flu_all_or, "flu_all_or_data.csv")

