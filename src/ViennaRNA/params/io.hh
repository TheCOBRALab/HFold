#ifndef IO
#define IO
#include <cstdio>


/**
 *  @brief  Default Energy Parameter File format
 *  @see    vrna_params_load(), vrna_params_load_from_string(),
 *          vrna_params_save()
 */
#define VRNA_PARAMETER_FORMAT_DEFAULT 0
/**
 *  @brief  Extract the filename from a file path
 */
char *vrna_basename(const char *path);

/**
 *  @brief Read a line of arbitrary length from a stream
 *
 *  Returns a pointer to the resulting string. The necessary memory is
 *  allocated and should be released using @e free() when the string is
 *  no longer needed.
 *
 *  @param  fp  A file pointer to the stream where the function should read from
 *  @return     A pointer to the resulting string
 */
char *vrna_read_line(FILE *fp);

/**
 *  @brief Load energy parameters from a file
 *
 *  @see  vrna_params_load_from_string(), vrna_params_save(),
 *        vrna_params_load_defaults(), vrna_params_load_RNA_Turner2004(),
 *        vrna_params_load_RNA_Turner1999(), vrna_params_load_RNA_Andronescu2007(),
 *        vrna_params_load_RNA_Langdon2018(), vrna_params_load_RNA_misc_special_hairpins(),
 *        vrna_params_load_DNA_Mathews2004(), vrna_params_load_DNA_Mathews1999()
 *
 *  @param  fname   The path to the file containing the energy parameters
 *  @param  options File format bit-mask (usually #VRNA_PARAMETER_FORMAT_DEFAULT)
 *  @return         Non-zero on success, 0 on failure
 */
int vrna_params_load(const char fname[], unsigned int options);

/**
 *  @brief  Load energy paramters from string
 *
 *  The string must follow the default energy parameter file convention!
 *  The optional @p name argument allows one to specify a name for the
 *  parameter set which is stored internally.
 *
 *  @see  vrna_params_load(), vrna_params_save(),
 *        vrna_params_load_defaults(), vrna_params_load_RNA_Turner2004(),
 *        vrna_params_load_RNA_Turner1999(), vrna_params_load_RNA_Andronescu2007(),
 *        vrna_params_load_RNA_Langdon2018(), vrna_params_load_RNA_misc_special_hairpins(),
 *        vrna_params_load_DNA_Mathews2004(), vrna_params_load_DNA_Mathews1999()
 *
 *
 *  @param    string  A 0-terminated string containing energy parameters
 *  @param    name    A name for the parameter set in @p string (Maybe @p NULL)
 *  @param    options File format bit-mask (usually #VRNA_PARAMETER_FORMAT_DEFAULT)
 *  @return           Non-zero on success, 0 on failure
 */
int vrna_params_load_from_string(const char *string, const char *name, unsigned int options);

/**
 *  @brief  Load default RNA energy parameter set
 *
 *  This is a convenience function to load the Turner 2004 RNA free
 *  energy parameters. It's the same as calling vrna_params_load_RNA_Turner2004()
 *
 *  @see  vrna_params_load(), vrna_params_load_from_string(),
 *        vrna_params_save(), vrna_params_load_RNA_Turner2004(),
 *        vrna_params_load_RNA_Turner1999(), vrna_params_load_RNA_Andronescu2007(),
 *        vrna_params_load_RNA_Langdon2018(), vrna_params_load_RNA_misc_special_hairpins(),
 *        vrna_params_load_DNA_Mathews2004(), vrna_params_load_DNA_Mathews1999()
 *
 *
 *  @return Non-zero on success, 0 on failure
 */
int vrna_params_load_defaults(void);

/**
 *  @brief  Load energy paramters from string
 *
 *  The string must follow the default energy parameter file convention!
 *  The optional @p name argument allows one to specify a name for the
 *  parameter set which is stored internally.
 *
 *  @see  vrna_params_load(), vrna_params_save(),
 *        vrna_params_load_defaults(), vrna_params_load_RNA_Turner2004(),
 *        vrna_params_load_RNA_Turner1999(), vrna_params_load_RNA_Andronescu2007(),
 *        vrna_params_load_RNA_Langdon2018(), vrna_params_load_RNA_misc_special_hairpins(),
 *        vrna_params_load_DNA_Mathews2004(), vrna_params_load_DNA_Mathews1999()
 *
 *
 *  @param    string  A 0-terminated string containing energy parameters
 *  @param    name    A name for the parameter set in @p string (Maybe @p NULL)
 *  @param    options File format bit-mask (usually #VRNA_PARAMETER_FORMAT_DEFAULT)
 *  @return           Non-zero on success, 0 on failure
 */
int vrna_params_load_from_string(const char *string, const char *name, unsigned int options);


/**
 *  @brief  Load Turner 2004 RNA energy parameter set
 *
 *  @see  vrna_params_load(), vrna_params_load_from_string(),
 *        vrna_params_save(), vrna_params_load_defaults(),
 *        vrna_params_load_RNA_Turner1999(), vrna_params_load_RNA_Andronescu2007(),
 *        vrna_params_load_RNA_Langdon2018(), vrna_params_load_RNA_misc_special_hairpins(),
 *        vrna_params_load_DNA_Mathews2004(), vrna_params_load_DNA_Mathews1999()
 *
 *
 *  @return Non-zero on success, 0 on failure
 */
int vrna_params_load_RNA_Turner2004(void);

/**
 *  @brief  Load Mathews 2004 DNA energy parameter set
 *
 *  @see  vrna_params_load(), vrna_params_load_from_string(),
 *        vrna_params_save(), vrna_params_load_RNA_Turner2004(),
 *        vrna_params_load_RNA_Turner1999(), vrna_params_load_RNA_Andronescu2007(),
 *        vrna_params_load_RNA_Langdon2018(), vrna_params_load_RNA_misc_special_hairpins(),
 *        vrna_params_load_defaults(), vrna_params_load_DNA_Mathews1999()
 *
 *
 *  @return Non-zero on success, 0 on failure
 */
int vrna_params_load_DNA_Mathews2004(void);

/**
 *  @brief Get the file name of the parameter file that was most recently loaded
 *
 *  @return The file name of the last parameter file, or NULL if parameters are still at defaults
 */
const char* last_parameter_file(void);

/**
 *  @brief
 *
 */
enum parset {
  UNKNOWN= -1, QUIT,
  S, S_H, HP, HP_H, B, B_H, IL, IL_H, MMH, MMH_H, MMI, MMI_H,
  MMI1N, MMI1N_H, MMI23, MMI23_H, MMM, MMM_H, MME, MME_H, D5, D5_H, D3, D3_H,
  INT11, INT11_H, INT21, INT21_H, INT22, INT22_H, ML, TL,
  TRI, HEX, NIN, MISC
};

/**
 *  @brief
 *
 */
enum parset gettype(const char *ident);


/**
 *  @brief
 *
 */
char const* settype(enum parset s);

#endif