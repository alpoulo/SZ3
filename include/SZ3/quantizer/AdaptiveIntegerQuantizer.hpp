#ifndef _SZ_ADAPTIVE_INTEGER_QUANTIZER_HPP
#define _SZ_ADAPTIVE_INTEGER_QUANTIZER_HPP

#include <cstring>
#include <cassert>
#include <iostream>
#include <vector>
#include "SZ3/def.hpp"
#include "SZ3/quantizer/Quantizer.hpp"

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
            adaptive_flag = (1 << ab);
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

        void set_max_error(double eb) {
            max_error = eb;
            for (int i = 0; i < adaptive_bits; i++) {
                max_error /= 2;
            }
        }

        double get_max_error() { return max_error; }

        std::vector<int> get_adapt_inds() { return adapt_inds; }

        // quantize the data with a prediction value, and returns the quantization index
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

        // quantize the data with a prediction value, and returns the quantization index and the decompressed data
        // int quantize(T data, T pred, T& dec_data);
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
                } else {
                    if (fabs(decompressed_data - data) < this->max_error) {
                        data = decompressed_data;
                        adapt_inds.push_back(adaptive_flag);
                    } else {
                        //data = decompressed_data;
                        //adapt_inds.push_back(adaptive_flag);
                        adaptive_quantize_encode(data, decompressed_data);
                    }
                    return quant_index_shifted;
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
            actual_adapt_ctr++;
            data = decompressed_data;
            adapt_inds.push_back(adaptive_index);
        }

        T adaptive_quantize_decode(T data, int quant_index, int adapt_index) {
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
            T data = pred + 2 * (quant_index - this->radius) * this->error_bound;
            if (adapt_ctr < adapt_inds.size()) {
                int adapt_index = adapt_inds[adapt_ctr++]; 
                if (adapt_index == adaptive_flag) {
                    return data;
                }
                return adaptive_quantize_decode(data, quant_index, adapt_index);
            }
            else {
                return data;
            }
        }

        T recover_unpred() {
            return unpred[index++];
        }

        size_t size_est() {
            return unpred.size() * sizeof(T) + sizeof(int); //include sizeof(int) to save # adaptive bits
            //return unpred.size() * sizeof(T);
        }

        void save(unsigned char *&c) const {
            // std::string serialized(sizeof(uint8_t) + sizeof(T) + sizeof(int),0);
            c[0] = 0b00000011;
            c += 1;
            *reinterpret_cast<int *>(c) = this->adaptive_bits;
            c += sizeof(int);
            *reinterpret_cast<double *>(c) = this->error_bound;
            c += sizeof(double);
            *reinterpret_cast<int *>(c) = this->radius;
            c += sizeof(int);
            *reinterpret_cast<size_t *>(c) = unpred.size();
            c += sizeof(size_t);
            memcpy(c, unpred.data(), unpred.size() * sizeof(T));
            c += unpred.size() * sizeof(T);
            //std::cout << "actual adaptive quantizations: " << actual_adapt_ctr << "\n";
        };

        void load(const unsigned char *&c, size_t &remaining_length) {
            assert(remaining_length > (sizeof(uint8_t) + sizeof(T) + sizeof(int)));
            c += sizeof(uint8_t);
            remaining_length -= sizeof(uint8_t);
            this->adaptive_bits = *reinterpret_cast<const int *>(c);
            this->adaptive_flag = 1 << adaptive_bits;
            c += sizeof(int);
            this->error_bound = *reinterpret_cast<const double *>(c);
            this->error_bound_reciprocal = 1.0 / this->error_bound;
            c += sizeof(double);
            this->radius = *reinterpret_cast<const int *>(c);
            c += sizeof(int);
            size_t unpred_size = *reinterpret_cast<const size_t *>(c);
            c += sizeof(size_t);
            this->unpred = std::vector<T>(reinterpret_cast<const T *>(c), reinterpret_cast<const T *>(c) + unpred_size);
            c += unpred_size * sizeof(T);
            // std::cout << "loading: eb = " << this->error_bound << ", unpred_num = "  << unpred.size() << std::endl;
            // reset index
            index = 0;
            adapt_ctr = 0;
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

        void postcompress_data(std::vector<int> &quant_inds) {
            quant_inds.insert(quant_inds.end(), adapt_inds.begin(), adapt_inds.end());
        }

        virtual void postdecompress_data() {
        }

        virtual void precompress_data() {};

        virtual void predecompress_data() {};

        void predecompress_data(std::vector<int> &quant_inds, size_t num_elements) {
            copy(quant_inds.begin() + num_elements, quant_inds.end(), back_inserter(adapt_inds));
            //std::cout << "num elements: " << num_elements << "\n";
            //std::cout << "adapt inds size: " << adapt_inds.size()  << "\n";
            //std::cout << "adapt_inds[0]: " << adapt_inds[0] << "\n";
        }


        int quantize_and_overwrite(T ori, T pred, T &dest) {
            T diff = ori - pred;
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
                if (fabs(decompressed_data - ori) > this->error_bound) {
                    unpred.push_back(ori);
                    dest = ori;
                    return 0;
                } else {
                    if (fabs(decompressed_data - ori) < max_error) {
                        dest = decompressed_data;
                        adapt_inds.push_back(adaptive_flag);
                        return quant_index_shifted;
                    } else {
                        dest = ori;
                        adaptive_quantize_encode(dest, decompressed_data);
                        return quant_index_shifted;
                    }
                }
            } else {
                unpred.push_back(ori);
                dest = ori;
                return 0;
            }
        }

    private:
        std::vector<T> unpred;
        size_t index = 0; // used in decompression only
        std::vector<int> adapt_inds;
        double error_bound;
        double error_bound_reciprocal;
        int radius; // quantization interval radius
        int adaptive_bits;
        int adaptive_flag;
        double max_error;
        size_t adapt_ctr = 0; // used in decompression only
        int actual_adapt_ctr = 0;
    };

}
#endif
