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
        AdaptiveLinearQuantizer() : error_bound(1), error_bound_reciprocal(1), radius(32768), adaptive_bits(2) {
        }

        AdaptiveLinearQuantizer(double eb, int ab, int r = 32768) : error_bound(eb),
                                                    error_bound_reciprocal(1.0 / eb),
                                                    radius(r), 
                                                    adaptive_bits(ab) {
            assert(eb != 0);
            set_max_error(eb);
        }

        int get_radius() const { return radius; }

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

        size_t get_num_adaptive() const { return actual_adapt_ctr; }

        std::vector<int> get_adapt_inds() { return adapt_inds; }

        /*
         * quantize data with a prediction value, return the quantization index and the decompressed data
         * uses value of '1' in quantization array to indicate adaptive quantization has occurred
         */
        int quantize_and_overwrite(T &data, T pred) {
            T diff = data - pred;
            int quant_index = (int) (fabs(diff) * this->error_bound_reciprocal) + 1;
            if (quant_index < this->radius * 2) {
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
                    unpred.push_back(data);
                    return 0;
                } else if (quant_index_shifted == 1) {
                    // value of '1' is reserved
                    unpred.push_back(data);
                    return 0;
                } else {
                    single_quant_ctr++;
                    if (fabs(decompressed_data - data) < this->max_error) {
                        data = decompressed_data;
                        return quant_index_shifted;
                    } else {
                        multi_quant_ctr++;
                        actual_inds.push_back(quant_index_shifted);
                        adaptive_quantize_encode(data, decompressed_data);
                        // save reserved value to indicate adaptive quantization occurred
                        return 1;
                    }
                }
            } else {
                unpred.push_back(data);
                return 0;
            }
        }

        void adaptive_quantize_encode(T &data, T decompressed_data) {
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
            //actual_adapt_ctr++;
            data = decompressed_data;
            actual_adapt_ctr++;
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
            if (quant_index == 1) {
                assert(actual_ctr < actual_inds.size());
                quant_index = actual_inds[actual_ctr++];
                T data = pred + 2 * (quant_index - this->radius) * this->error_bound;
                int adapt_index = adapt_inds[adapt_ctr++];
                return adaptive_quantize_decode(data, adapt_index);
            }
            else {
                T data = pred + 2 * (quant_index - this->radius) * this->error_bound;
                return data;
            }
        }

        T recover_unpred() {
            return unpred[index++];
        }

        size_t size_est() {
            //include sizeof(int) to save # adaptive bits
            return unpred.size() * sizeof(T) + sizeof(int) * 3; 
        }

        void save(unsigned char *&c) const {
            uchar *buffer_pos = c;

            // std::string serialized(sizeof(uint8_t) + sizeof(T) + sizeof(int),0);
            c[0] = 0b00000011;
            c += 1;
            *reinterpret_cast<double *>(c) = this->error_bound;
            c += sizeof(double);
            *reinterpret_cast<int *>(c) = this->radius;
            c += sizeof(int);
            *reinterpret_cast<int *>(c) = this->adaptive_bits;
            c += sizeof(int);
            *reinterpret_cast<int *>(c) = this->actual_adapt_ctr;
            c += sizeof(int);
            *reinterpret_cast<size_t *>(c) = unpred.size();
            c += sizeof(size_t);
            memcpy(c, unpred.data(), unpred.size() * sizeof(T));
            c += unpred.size() * sizeof(T);

            size_t outsize = 0;

            if (!adapt_inds.empty()) {
                HuffmanEncoder<int> encoder = HuffmanEncoder<int>();
                encoder.preprocess_encode(adapt_inds, 0);
                encoder.save(c);
                outsize = encoder.encode(adapt_inds, c);
                encoder.postprocess_encode();
            }

            std::cout << "outsize = " << outsize << std::endl;
            std::cout << "pred values = " << single_quant_ctr << std::endl;
            std::cout << "adapt values = " << multi_quant_ctr << std::endl;
        };

        void load(const unsigned char *&c, size_t &remaining_length) {
            //assert(remaining_length > (sizeof(uint8_t) + sizeof(T) + sizeof(int)));
            const unsigned char *buffer_pos = c;

            c += sizeof(uint8_t);
            remaining_length -= sizeof(uint8_t);
            this->error_bound = *reinterpret_cast<const double *>(c);
            this->error_bound_reciprocal = 1.0 / this->error_bound;
            c += sizeof(double);
            this->radius = *reinterpret_cast<const int *>(c);
            c += sizeof(int);
            this->adaptive_bits = *reinterpret_cast<const int *>(c);
            c += sizeof(int);
            this->actual_adapt_ctr = *reinterpret_cast<const int *>(c);
            c += sizeof(int);
            size_t unpred_size = *reinterpret_cast<const size_t *>(c);
            c += sizeof(size_t);
            this->unpred = std::vector<T>(reinterpret_cast<const T *>(c), reinterpret_cast<const T *>(c) + unpred_size);
            c += unpred_size * sizeof(T);

            if(actual_adapt_ctr != 0) {
                HuffmanEncoder<int> decoder = HuffmanEncoder<int>();
                decoder.load(c, remaining_length);
                auto adapt_inds = decoder.decode(c, actual_adapt_ctr);
                this->adapt_inds = adapt_inds;
                decoder.postprocess_decode();
                remaining_length -= actual_adapt_ctr * sizeof(int);
            }
            
            // std::cout << "loading: eb = " << this->error_bound << ", unpred_num = "  << unpred.size() << std::endl;
            // reset index
            index = 0;
            adapt_ctr = 0;
            actual_ctr = 0;
        }

        void print() {
            printf("[AdaptiveIntegerQuantizer] error_bound = %.8G, radius = %d, unpred = %lu\n", error_bound, radius, unpred.size());
        }

        void clear() {
            unpred.clear();
            index = 0;
        }

        virtual void postcompress_data() {
        }

        /*
         * pack adaptive bits to save space and losslessly store with unpredictable data
         * use value of '1' in quant_inds to indicate adaptive quantization occurred
         * return true quantization indices (appended) with quant_inds to be huffman encoded
         * quant_inds.size() = num_elements + # of elts adaptively quantized
         */
        void postcompress_data(std::vector<int> &quant_inds) {
            quant_inds.insert(quant_inds.end(), actual_inds.begin(), actual_inds.end());

            char path[BUFSIZE];
            const char *envvar = "ADAPTENTROPY";
            if (!getenv(envvar)) {
                fprintf(stderr, "the environment variable %s was not found!\n", envvar);
                exit(1);
            }
            if (snprintf(path, BUFSIZE, "%s", getenv(envvar)) >= BUFSIZE) {
                fprintf(stderr, "BUFSIZE of %d was too small! aborting.\n", BUFSIZE);
                exit(1);
            }
            writefile(path, adapt_inds.data(), adapt_inds.size());
        }
        
        virtual void postdecompress_data() {
        }

        virtual void precompress_data() {};

        virtual void predecompress_data() {};

        /*
         * get true quantization indices of adaptively quantized data from quantization index array
         * unpack adaptive quantization bits
         * quant_inds.size() should be num_elements + # elts adaptively quantized when passed in
         */
        void predecompress_data(std::vector<int> &quant_inds, size_t num_elements) {
            copy(quant_inds.begin() + num_elements, quant_inds.end(), back_inserter(actual_inds));
        }

        int quantize(T data, T pred) {
            T diff = data - pred;
            int quant_index = (int) (fabs(diff) * this->error_bound_reciprocal) + 1;
            if (quant_index < this->radius * 2) {
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
                    return 0;
                } else {
                    return quant_index_shifted;
                }
            } else {
                return 0;
            }
        }

    private:
        std::vector<T> unpred;
        size_t index = 0; // used in decompression only
        std::vector<int> adapt_inds;
        std::vector<int> actual_inds;
        //int adapt_ctr;
        double error_bound;
        double error_bound_reciprocal;
        int radius; // quantization interval radius
        int adaptive_bits;
        double max_error;
        size_t adapt_ctr = 0; // used in decompression only
        int actual_adapt_ctr = 0;
        int packed_size = 0;
        int actual_ctr = 0;
        int single_quant_ctr = 0;
        int multi_quant_ctr = 0;
    };

}
#endif
