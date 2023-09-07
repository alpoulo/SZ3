#ifndef _SZ_ADAPTIVE_INTEGER_QUANTIZER_HPP
#define _SZ_ADAPTIVE_INTEGER_QUANTIZER_HPP

#include <cstring>
#include <cassert>
#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include "SZ3/def.hpp"
#include "SZ3/quantizer/Quantizer.hpp"
#include "SZ3/encoder/HuffmanEncoder.hpp"
#include "SZ3/utils/FileUtil.hpp"

#define BUFSIZE 256
namespace SZ {

    template<class T>
    class AdaptiveLinearQuantizer : public concepts::QuantizerInterface<T> {
    public:
        AdaptiveLinearQuantizer() : error_bound(1), error_bound_reciprocal(1), radius(32768), adaptive_bits(2), aqMode(RA_HUFF) {
        }

        AdaptiveLinearQuantizer(double eb, int ab, uint8_t mode, int r = 32768) : 
                                                    error_bound(eb),
                                                    error_bound_reciprocal(1.0 / eb),
                                                    radius(r), 
                                                    adaptive_bits(ab), aqMode(mode) {
          assert(eb != 0);
          set_max_error(eb);
          adaptive_flag = 1 << ab;
        }

        AdaptiveLinearQuantizer(double eb, int ab, uint8_t mode, bool hist_flag, int r = 32768) : 
                                                    error_bound(eb),
                                                    error_bound_reciprocal(1.0 / eb),
                                                    radius(r), 
                                                    adaptive_bits(ab), aqMode(mode) {
          assert(eb != 0);
          this->hist_flag = hist_flag;
          set_max_error(eb);
          adaptive_flag = 1 << ab;
          if (hist_flag) {
            if (aqMode == RA_APPEND || aqMode == RA_PACK || aqMode == RA_HUFF) {
              adapthist_size = adaptive_flag + 1;
            }
            else {
              adapthist_size = adaptive_flag;
            }
            adapthist = (int *)malloc(sizeof(int) * adapthist_size);
            memset(adapthist,0, sizeof (int) * adapthist_size);
            hist_size = 2 * radius;
            this->hist.resize(hist_size);
            std::fill(this->hist.begin(), this->hist.end(), 0);
          }
        }

        void update_predictability() {
          this->pred_rate = (double)this->single_quant_ctr / (double)this->total_points;
          if (!use_adapt) {
            if (pred_rate < adapt_thres) {
              use_adapt = true;
            }
          } else {
            if (pred_rate >= adapt_thres) {
              use_adapt = false;
            }
          }
        }

        int get_radius() const { return radius; }

        size_t get_num_elements(size_t num_elements) const { 
          if (aqMode == RS_PACK || aqMode == RS_HUFF) {
            return num_elements + multi_quant_ctr;
          }
          else if (aqMode == RA_APPEND) {
            return num_elements + single_quant_ctr;
          } 
          else {
            return num_elements;
          }
        }
        
        double get_eb() const { return error_bound; }

        void set_eb(double eb) {
          error_bound = eb;
          error_bound_reciprocal = 1.0 / eb;
        }

        int get_ab() const { return adaptive_bits; }

        void set_ab(int ab) {
          adaptive_bits = ab;
        }

        double get_max_error() { return max_error; }

        void set_max_error(double eb) {
          max_error = eb;
          for (int i = 0; i < adaptive_bits; i++) {
              max_error /= 2;
          }
        }

        size_t get_num_adaptive() const { return multi_quant_ctr; }

        std::vector<int> get_adapt_inds() { return adapt_inds; }

        void update_hist(int quant_index) {
          if (hist_flag) {
              this->hist[quant_index]++;
          }
        }

        void update_adapthist(int quant_index) {
          if (hist_flag) {
              this->adapthist[quant_index]++;
          }
        }

        int quantize(T data, T pred) {
          T diff = data - pred;
          int quant_index = (int) (fabs(diff) * this->error_bound_reciprocal) + 1;
          //if (quant_index < this->radius * 2) {
          if ((quant_index < this->radius * 2) && (quant_index >= 0)) {
            quant_index >>= 1;
            int half_index = quant_index;
            quant_index <<= 1;
            int quant_index_shifted;
            if (diff < 0) {
              quant_index = -quant_index;
              quant_index_shifted = this->radius - half_index;
            } else {
                quant_index_shifted = this->radius + half_index;
            }
            T decompressed_data = pred + quant_index * this->error_bound;
            if (fabs(decompressed_data - data) > this->error_bound) {
              // unpredictable data
              unpred.push_back(data);  
              update_hist(0);
              return 0;
            }
            if ((aqMode == RS_PACK || aqMode == RS_HUFF) && quant_index_shifted == 1) {
              unpred.push_back(data);
              update_hist(0);
              return 0;
            }
            single_quant_ctr++;
            update_hist(quant_index_shifted);
           
            // not worth it to use AQ
            if (fabs(decompressed_data - data) < this->max_error) {
              if (aqMode == RA_APPEND || aqMode == RA_PACK || aqMode == RA_HUFF) {  
                // adaptive flag indicates no adaptive quantization occurred
                adapt_inds.push_back(adaptive_flag);
              }
            } else {
              multi_quant_ctr++;
              adaptive_quantize_encode(data, decompressed_data);
              
              if (aqMode == RS_PACK || aqMode == RS_HUFF) {
                actual_inds.push_back(quant_index_shifted);
                update_hist(1);
                return 1;
              }
            }
            return quant_index_shifted;
          } 
          else {
            unpred.push_back(data);
            update_hist(0);
            update_predictability();
            return 0;
          }
        }

        /*
         * quantize data with a prediction value, return the quantization index and the decompressed data
         * uses value of '1' in quantization array to indicate adaptive quantization has occurred
         */
        int quantize_and_overwrite(T &data, T pred) {
          T diff = data - pred;
          int quant_index = (int) (fabs(diff) * this->error_bound_reciprocal) + 1;
          //if (quant_index < this->radius * 2) {
          if ((quant_index < this->radius * 2) && (quant_index >= 0)) {
            quant_index >>= 1;
            int half_index = quant_index;
            quant_index <<= 1;
            int quant_index_shifted;
            if (diff < 0) {
              quant_index = -quant_index;
              quant_index_shifted = this->radius - half_index;
            } else {
              quant_index_shifted = this->radius + half_index;
            }
            T decompressed_data = pred + quant_index * this->error_bound;
            total_points++;
            if (fabs(decompressed_data - data) > this->error_bound) {
              unpred.push_back(data);
              update_hist(0);
              return 0;
            }
            if ((aqMode == RS_PACK || aqMode == RS_HUFF) && quant_index_shifted == 1) {
              unpred.push_back(data);
              update_hist(0);
              return 0;
            }
            single_quant_ctr++;
            update_hist(quant_index_shifted);
            
            // not worth it to use AQ
            if (fabs(decompressed_data - data) < this->max_error) {
              data = decompressed_data;
              if (aqMode == RA_APPEND || aqMode == RA_PACK || aqMode == RA_HUFF) {
                // adaptive flag indicates no AQ occurred
                adapt_inds.push_back(adaptive_flag);
                update_adapthist(adaptive_flag);
              }
            } else {
              multi_quant_ctr++;
              adaptive_quantize_encode_and_overwrite(data, decompressed_data);
              
              if (aqMode == RS_PACK || aqMode ==RS_HUFF) {
                actual_inds.push_back(quant_index_shifted);
                update_hist(1);
                return 1;
              }
            }
            return quant_index_shifted;
          } else {
            unpred.push_back(data);
            update_hist(0);
            return 0;
          }
        }

        /*
         * quantize data with a prediction value, return the quantization index and the decompressed data
         * uses value of '1' in quantization array to indicate adaptive quantization has occurred
         */
        int quantize_and_overwrite(T ori, T pred, T &dest) {
          T diff = ori - pred;
          int quant_index = (int) (fabs(diff) * this->error_bound_reciprocal) + 1;
          if ((quant_index < this->radius * 2) && (quant_index >= 0)) {
            quant_index >>= 1;
            int half_index = quant_index;
            quant_index <<= 1;
            int quant_index_shifted;
            if (diff < 0) {
              quant_index = -quant_index;
              quant_index_shifted = this->radius - half_index;
            } else {
              quant_index_shifted = this->radius + half_index;
            }
            T decompressed_data = pred + quant_index * this->error_bound;
            if (fabs(decompressed_data - ori) > this->error_bound) {
              unpred.push_back(ori);
              dest = ori;
              update_hist(0);
              return 0;
            } 
            if ((aqMode == RS_PACK || aqMode == RS_HUFF) && quant_index_shifted == 1) {
              unpred.push_back(ori);
              dest = ori;
              update_hist(0);
              return 0;
            }
            single_quant_ctr++;
            update_hist(quant_index_shifted);

            if (fabs(decompressed_data - ori) < this->max_error) {
              if (aqMode == RA_APPEND || aqMode == RA_PACK || aqMode == RA_HUFF) {
                // adaptive flag indicates no adaptive quantization occurred
                adapt_inds.push_back(adaptive_flag);
              }
              dest = decompressed_data;
              update_adapthist(adaptive_flag);
            } else {
              multi_quant_ctr++;
              adaptive_quantize_encode_and_overwrite(ori, decompressed_data, dest);

              if (aqMode == RS_PACK || aqMode == RS_HUFF) {
                actual_inds.push_back(quant_index_shifted);
                update_hist(1);
                return 1;
              }
            }
            return quant_index_shifted;
          } else {
            unpred.push_back(ori);
            dest = ori;
            update_hist(0);
            return 0;
          }
        }

        void adaptive_quantize_encode(T data, T decompressed_data) {
          double precision = this->error_bound;
          int adaptive_index = 0;
          for (int i = 0; i < adaptive_bits; i++) {
            precision /= 2;
            if ((data - decompressed_data) > 0) {
              adaptive_index ^= (1 << i);
              decompressed_data += precision;
            } else {
              decompressed_data -= precision;
            }
          }
          update_adapthist(adaptive_index);
          adapt_inds.push_back(adaptive_index);
        }

        void adaptive_quantize_encode_and_overwrite(T &data, T decompressed_data) {
          double precision = this->error_bound;
          int adaptive_index = 0;
          for (int i = 0; i < adaptive_bits; i++) {
            precision /= 2;
            if ((data - decompressed_data) > 0) {
              adaptive_index ^= (1 << i);
              decompressed_data += precision;
            } else {
              decompressed_data -= precision;
            }
          }
          data = decompressed_data;
          adapt_inds.push_back(adaptive_index);
          update_adapthist(adaptive_index);
        }

        void adaptive_quantize_encode_and_overwrite(T ori, T decompressed_data, T &dest) {
          double precision = this->error_bound;
          int adaptive_index = 0;
          for (int i = 0; i < adaptive_bits; i++) {
            precision /= 2;
            if ((ori - decompressed_data) > 0) {
              adaptive_index ^= (1 << i);
              decompressed_data += precision;
            } else {
              decompressed_data -= precision;
            }
          }
          dest = decompressed_data;
          update_adapthist(adaptive_index);
          adapt_inds.push_back(adaptive_index);
        }

        T adaptive_quantize_decode(T data, int adapt_index) {
          T decompressed_data = data;
          double precision = error_bound;
          for (int i = 0; i < adaptive_bits; i++) {
            precision /= 2;
            if ((adapt_index >> i) & 1) {
              decompressed_data += precision;
            } else {
              decompressed_data -= precision;
            }
          }
          return decompressed_data;
        }

        // recover the data using the quantization index
        T recover(T pred, int quant_index) {
          if (quant_index) {
            return recover_pred(pred, quant_index);
          } else {
            return recover_unpred();
          }
        }

        T recover_pred(T pred, int quant_index) {
          if ((aqMode == RS_PACK || aqMode == RS_HUFF) && quant_index == 1) {
            assert(actual_ctr < actual_inds.size());
            quant_index = actual_inds[actual_ctr++];
            T data = pred + 2 * (quant_index - this->radius) * this->error_bound;
            if (adapt_ctr < adapt_inds.size()) {
              int adapt_index = adapt_inds[adapt_ctr++];
              return adaptive_quantize_decode(data, adapt_index);
            } else {
              return data;
            }
          }
          else {
            T data = pred + 2 * (quant_index - this->radius) * this->error_bound;
            if (adapt_ctr < adapt_inds.size()) {
              int adapt_index = adapt_inds[adapt_ctr++];
              if (adapt_index == adaptive_flag) {
                return data;
              }
              return adaptive_quantize_decode(data, adapt_index);
            } else {
              return data;
            }
          }
        }

        T recover_unpred() {
          return unpred[index++];
        }

        size_t size_est() {
          //include sizeof(int) to save # adaptive bits
          if (aqMode == RS_PACK || aqMode == RA_PACK) {
            return unpred.size() * sizeof(T) + packed_inds.size() * sizeof(int) + sizeof(uint8_t) + sizeof(int) * 4;
          }
          return unpred.size() * sizeof(T) + sizeof(uint8_t) + sizeof(int) * 4; 
        }

        void save(unsigned char *&c) const {
          // std::string serialized(sizeof(uint8_t) + sizeof(T) + sizeof(int),0);
          c[0] = 0b00000011;
          c += 1;
          *reinterpret_cast<double *>(c) = this->error_bound;
          c += sizeof(double);
          *reinterpret_cast<int *>(c) = this->radius;
          c += sizeof(int);
          *reinterpret_cast<int *>(c) = this->adaptive_bits; // # adaptive bits
          c += sizeof(int);
          *reinterpret_cast<uint8_t *>(c) = this->aqMode;
          c += sizeof(uint8_t);
          *reinterpret_cast<int *>(c) = single_quant_ctr;
          c += sizeof(int);
          *reinterpret_cast<int *>(c) = multi_quant_ctr;
          c += sizeof(int);
          *reinterpret_cast<size_t *>(c) = unpred.size();
          c += sizeof(size_t);
          memcpy(c, unpred.data(), unpred.size() * sizeof(T));
          c += unpred.size() * sizeof(T);

          if (aqMode == RS_PACK || aqMode == RA_PACK) {
            *reinterpret_cast<int *>(c) = this->packed_inds.size();
            c += sizeof(int);
            memcpy(c, packed_inds.data(), packed_inds.size() * sizeof(int));
            c += packed_inds.size() * sizeof(int);
            std::cout << "quantPackedOutsize = " << packed_inds.size() << std::endl;
          }
          
          if (aqMode == RS_HUFF || aqMode == RA_HUFF) {
            size_t outsize = 0;
            if (!adapt_inds.empty()) {
              //HuffmanEncoder<int> encoder = HuffmanEncoder<int>();
              HuffmanEncoder<int> encoder = HuffmanEncoder<int>();
              encoder.preprocess_encode(adapt_inds, 0);
              encoder.save(c);
              outsize = encoder.encode(adapt_inds, c);
              encoder.postprocess_encode();
            }
            std::cout << "quantHuffmanOutsize = " << outsize << std::endl;
          }
          if (aqMode == RA_APPEND) {
            std::cout << "quantHuffmanOutsize = 0" << std::endl;
          }
          
          std::cout << "pred values = " << single_quant_ctr << std::endl;
          std::cout << "adapt values = " << multi_quant_ctr << std::endl;
        };

        void save(unsigned char *&c, int indicator) const {
          // std::string serialized(sizeof(uint8_t) + sizeof(T) + sizeof(int),0);
          c[0] = 0b00000011;
          c += 1;
          *reinterpret_cast<double *>(c) = this->error_bound;
          c += sizeof(double);
          *reinterpret_cast<int *>(c) = this->radius;
          c += sizeof(int);
          *reinterpret_cast<int *>(c) = this->adaptive_bits;
          c += sizeof(int);
          *reinterpret_cast<uint8_t *>(c) = this->aqMode;
          c += sizeof(uint8_t);
          *reinterpret_cast<int *>(c) = this->single_quant_ctr;
          c += sizeof(int);
          *reinterpret_cast<int *>(c) = this->multi_quant_ctr;
          c += sizeof(int);
          *reinterpret_cast<size_t *>(c) = unpred.size();
          c += sizeof(size_t);
          memcpy(c, unpred.data(), unpred.size() * sizeof(T));
          c += unpred.size() * sizeof(T);

          print_bins(indicator);
          
          if (aqMode == RS_PACK || aqMode == RA_PACK) {
            *reinterpret_cast<int *>(c) = this->packed_inds.size();
            c += sizeof(int);
            memcpy(c, packed_inds.data(), packed_inds.size() * sizeof(int));
            c += packed_inds.size() * sizeof(int);
            std::cout << "quantPackedOutsize = " << packed_inds.size() * sizeof(int) << std::endl;
          }
          if (aqMode == RS_HUFF || aqMode == RA_HUFF) {
            size_t outsize = 0;
            if (!adapt_inds.empty()) {
              //HuffmanEncoder<int> encoder = HuffmanEncoder<int>();
              HuffmanEncoder<int> encoder = HuffmanEncoder<int>();
              encoder.preprocess_encode(adapt_inds, 0);
              encoder.save(c);
              outsize = encoder.encode(adapt_inds, c);
              encoder.postprocess_encode();
            }
            std::cout << "quantHuffmanOutsize = " << outsize << std::endl;
          }
          if (aqMode == RA_APPEND) {
            std::cout << "quantHuffmanOutsize = 0" << std::endl;
          }
          std::cout << "pred values = " << single_quant_ctr << std::endl;
          std::cout << "adapt values = " << multi_quant_ctr << std::endl;
        };

        void load(const unsigned char *&c, size_t &remaining_length) {
          //assert(remaining_length > (sizeof(uint8_t) + sizeof(T) + sizeof(int)));
          c += sizeof(uint8_t);
          remaining_length -= sizeof(uint8_t);
          this->error_bound = *reinterpret_cast<const double *>(c);
          this->error_bound_reciprocal = 1.0 / this->error_bound;
          c += sizeof(double);
          this->radius = *reinterpret_cast<const int *>(c);
          c += sizeof(int);
          this->adaptive_bits = *reinterpret_cast<const int *>(c);
          c += sizeof(int);
          this->aqMode = *reinterpret_cast<const uint8_t *>(c);
          c += sizeof(uint8_t);
          this->adaptive_flag = (1 << adaptive_bits);
          this->single_quant_ctr = *reinterpret_cast<const int *>(c);
          c += sizeof(int);
          this->multi_quant_ctr = *reinterpret_cast<const int *>(c);
          c += sizeof(int);
          size_t unpred_size = *reinterpret_cast<const size_t *>(c);
          c += sizeof(size_t);
          this->unpred = std::vector<T>(reinterpret_cast<const T *>(c), reinterpret_cast<const T *>(c) + unpred_size);
          c += unpred_size * sizeof(T);
          
          if (aqMode == RS_PACK || aqMode == RA_PACK) {
            this->packed_size = *reinterpret_cast<const int *>(c);
            c += sizeof(int);
            this->packed_inds = std::vector<int>(reinterpret_cast<const int *>(c), reinterpret_cast<const int *>(c) + packed_size);
            c += packed_size * sizeof(int);
            if (aqMode == RA_PACK) {
              adaptive_unpack();
            }
          }

          if (aqMode == RS_HUFF || aqMode == RA_HUFF) {
            if(multi_quant_ctr != 0) {
              //HuffmanEncoder<int> decoder = HuffmanEncoder<int>();
              HuffmanEncoder<int> decoder = HuffmanEncoder<int>();
              decoder.load(c, remaining_length);
              auto adapt_inds = decoder.decode(c, multi_quant_ctr);
              this->adapt_inds = adapt_inds;
              decoder.postprocess_decode();
              remaining_length -= multi_quant_ctr * sizeof(int);
            }
          }
          // reset index
          index = 0;
          adapt_ctr = 0;
          actual_ctr = 0;
        }

        void print() {
          printf("[AdaptiveIntegerQuantizer] error_bound = %.8G, radius = %d, unpred = %lu\n", error_bound, radius, unpred.size());
        }

        void print_hist() const {
          char path[BUFSIZE];
          const char *envvar = "QUANTHISTOGRAM";
          if (!getenv(envvar)) {
            fprintf(stderr, "the environment variable %s was not found!\n", envvar);
            return;
          }
          if (snprintf(path, BUFSIZE, "%s", getenv(envvar)) >= BUFSIZE) {
            fprintf(stderr, "BUFSIZE of %d was too small! aborting.\n", BUFSIZE);
            return;
          }
          FILE *f = fopen(path, "w");
          if (f == NULL) {
            fprintf(stderr, "unable to open FILE %s for writing!\n", path);
            return;
          }
          int i;
          for (int i = 0; i < this->hist_size; i++) {
            fprintf(f, "%d ", hist[i]);
          }
          //writefile(path, adapt_inds.data(), adapt_inds.size());
          fclose(f);
          //free(hist);
        }

        void print_adapthist() const {
          char path[BUFSIZE];
          const char *envvar = "ADAPTHISTOGRAM";
          if (!getenv(envvar)) {
              fprintf(stderr, "the environment variable %s was not found!\n", envvar);
              return;
          }
          if (snprintf(path, BUFSIZE, "%s", getenv(envvar)) >= BUFSIZE) {
              fprintf(stderr, "BUFSIZE of %d was too small! aborting.\n", BUFSIZE);
              return;
          }
          FILE *f = fopen(path, "w");
          if (f == NULL) {
              fprintf(stderr, "unable to open FILE %s for writing!\n", path);
              return;
          }
          //fprintf(stderr, "hist_size: %d\n", hist_size);
          //fprintf(stderr, "fpath: %s\n", path);
          int i;
          for (int i = 0; i < this->adapthist_size; i++) {
              fprintf(f, "%d ", adapthist[i]);
          }
          //writefile(path, adapt_inds.data(), adapt_inds.size());
          fclose(f);
          free(adapthist);
        }

        void print_bins() const {
          char path[BUFSIZE];
          const char *envvar = "ADAPTENTROPY";
          if (!getenv(envvar)) {
              fprintf(stderr, "the environment variable %s was not found!\n", envvar);
              return;
          }
          if (snprintf(path, BUFSIZE, "%s", getenv(envvar)) >= BUFSIZE) {
              fprintf(stderr, "BUFSIZE of %d was too small! aborting.\n", BUFSIZE);
              return;
          }
          writefile(path, adapt_inds.data(), adapt_inds.size());
        }

        void print_bins(int indicator) const {
          char path[BUFSIZE];
          if (indicator == 0) {
              const char *envvar = "ADAPTENTROPYINDEPENDENT";
              if (!getenv(envvar)) {
                  fprintf(stderr, "the environment variable %s was not found!\n", envvar);
                  return;
              }
              if (snprintf(path, BUFSIZE, "%s", getenv(envvar)) >= BUFSIZE) {
                  fprintf(stderr, "BUFSIZE of %d was too small! aborting.\n", BUFSIZE);
                  return;
              }
              writefile(path, adapt_inds.data(), adapt_inds.size());
          } else { 
              const char *envvar = "ADAPTENTROPYLINER";
              if (!getenv(envvar)) {
                  fprintf(stderr, "the environment variable %s was not found!\n", envvar);
                  return;
              }
              if (snprintf(path, BUFSIZE, "%s", getenv(envvar)) >= BUFSIZE) {
                  fprintf(stderr, "BUFSIZE of %d was too small! aborting.\n", BUFSIZE);
                  return;
              }
              writefile(path, adapt_inds.data(), adapt_inds.size());
          } 
        }

        void adaptive_pack() {
          int shift = adaptive_bits;
          if (aqMode == RA_PACK) {
            shift = adaptive_bits + 1;
          }
          int max_shift = (sizeof(int) * 8) - shift;
          int offset = max_shift;
          int pack_value = 0;
          for (int i = 0; i < adapt_inds.size(); i++) {
            pack_value ^= (adapt_inds[i] << offset);
            offset -= shift;
            if(offset < 0) {
              packed_inds.push_back(pack_value);
              packed_size++;
              offset = max_shift;
              pack_value = 0;
            }
          }
          if(offset != max_shift) {
            packed_inds.push_back(pack_value);
            packed_size++;
          }
          std::cout << "packed_size: " << packed_size << std::endl;
        }

        void adaptive_unpack() {
          std::cout << "unpacking adaptive bits . . . ";
          int shift = adaptive_bits;
          if (aqMode == RA_PACK) {
            shift = adaptive_bits + 1;
          }
          int max_shift = (sizeof(int) * 8) - shift;
          int bit_check = 0;
          for (int i = 0; i < shift; i++) {
            bit_check ^= 1 << i;
          }
          for (int i = 0; i < packed_inds.size(); i++) {
            int cur_val = packed_inds[i];
            int offset = max_shift;
            while(offset >= 0) {
              adapt_inds.push_back((cur_val >> offset) & bit_check);
              offset -= shift;
            }
          }
          std::cout << "done\n";
        }

        void clear() {
          unpred.clear();
          adapt_inds.clear();
          hist.clear();
          index = 0;
        }

        void postcompress_data() {
          if (hist_flag) {
              print_hist();
              print_adapthist();
              print_bins();
          }
        }

        void postcompress_data(std::vector<int> &quant_inds) {
          if (hist_flag) {
              print_hist();
              print_adapthist();
              print_bins();
          }
          
          if (aqMode == RS_PACK || aqMode == RA_PACK) {
            adaptive_pack();
          }
          if (aqMode == RS_PACK || aqMode == RS_HUFF) {
            quant_inds.insert(quant_inds.end(), actual_inds.begin(), actual_inds.end());
          }
          if (aqMode == RA_APPEND) {
            quant_inds.insert(quant_inds.end(), adapt_inds.begin(), adapt_inds.end());
          }
        }
        
        virtual void postdecompress_data() {
        }

        virtual void precompress_data() {};

        virtual void predecompress_data() {};

        void predecompress_data(std::vector<int> &quant_inds, size_t num_elements) {
          if (aqMode == RS_PACK || aqMode == RS_HUFF) {
            copy(quant_inds.begin() + num_elements, quant_inds.end(), back_inserter(actual_inds));
          }
          if (aqMode == RS_PACK) {
            adaptive_unpack();
          }
          if (aqMode == RA_APPEND) {
            copy(quant_inds.begin() + num_elements, quant_inds.end(), back_inserter(adapt_inds));
          }
        }

    private:
      std::vector<T> unpred;
      size_t index = 0; // used in decompression only
      std::vector<int> adapt_inds;
      std::vector<int> packed_inds;
      std::vector<int> actual_inds;
      double error_bound;
      double error_bound_reciprocal;
      int radius; // quantization interval radius
      int adaptive_bits;
      int adaptive_flag = 0;
      uint8_t aqMode;
      double max_error;
      size_t adapt_ctr = 0; // used in decompression only
      int packed_size = 0;
      int actual_ctr = 0;
      int single_quant_ctr = 0;
      int multi_quant_ctr = 0;
      //int *hist;
      std::vector<int> hist;
      int hist_size;
      int hist_flag;
      int *adapthist;
      int adapthist_size;
      size_t true_num_elements;
      double pred_rate = 0.0;
      int total_points = 0;
      bool use_adapt = 0;
      double adapt_thres = 0.85;
      int most_used_bin = 0;
      int most_used_bin_amount = 0;
    };

}
#endif
