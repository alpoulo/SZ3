#ifndef SZ_ADAPTIVE_COMPRESSOR_HPP
#define SZ_ADAPTIVE_COMPRESSOR_HPP

#include "SZ3/compressor/Compressor.hpp"
#include "SZ3/frontend/Frontend.hpp"
#include "SZ3/encoder/Encoder.hpp"
#include "SZ3/lossless/Lossless.hpp"
#include "SZ3/utils/FileUtil.hpp"
#include "SZ3/utils/Config.hpp"
#include "SZ3/utils/Timer.hpp"
#include "SZ3/def.hpp"
#include <cstring>

namespace SZ {
    template<class T, uint N, class Frontend, class Encoder, class Lossless>
    class SZAdaptiveCompressor : public concepts::CompressorInterface<T> {
    public:


        SZAdaptiveCompressor(Frontend frontend, Encoder encoder, Lossless lossless) :
                frontend(frontend), encoder(encoder), lossless(lossless) {
            static_assert(std::is_base_of<concepts::FrontendInterface<T, N>, Frontend>::value,
                          "must implement the frontend interface");
            static_assert(std::is_base_of<concepts::EncoderInterface<int>, Encoder>::value,
                          "must implement the encoder interface");
            static_assert(std::is_base_of<concepts::LosslessInterface, Lossless>::value,
                          "must implement the lossless interface");
        }

        uchar *compress(const Config &conf, T *data, size_t &compressed_size) {


            std::vector<int> quant_inds = frontend.compress(data);
           
            //std::vector<int> adapt_inds(quant_inds.begin() + frontend.get_num_elements(), quant_inds.end());
            //quant_inds.resize(frontend.get_num_elements());

            encoder.preprocess_encode(quant_inds, 0);
            size_t bufferSize = 1.2 * (frontend.size_est() + encoder.size_est() + sizeof(T) * quant_inds.size());

            uchar *buffer = new uchar[bufferSize];
            uchar *buffer_pos = buffer;

            frontend.save(buffer_pos);

            encoder.save(buffer_pos);
            encoder.encode(quant_inds, buffer_pos);
            encoder.postprocess_encode();

            assert(buffer_pos - buffer < bufferSize);

            uchar *lossless_data = lossless.compress(buffer, buffer_pos - buffer, compressed_size);
            lossless.postcompress_data(buffer);

            return lossless_data;
        }

        T *decompress(uchar const *cmpData, const size_t &cmpSize, size_t num) {
            T *dec_data = new T[num];
            return decompress(cmpData, cmpSize, dec_data);
        }

        T *decompress(uchar const *cmpData, const size_t &cmpSize, T *decData) {
            size_t remaining_length = cmpSize;

            Timer timer(true);
            auto compressed_data = lossless.decompress(cmpData, remaining_length);
            uchar const *compressed_data_pos = compressed_data;
//            timer.stop("Lossless");

            frontend.load(compressed_data_pos, remaining_length);

            encoder.load(compressed_data_pos, remaining_length);

            timer.start();
            auto quant_inds = encoder.decode(compressed_data_pos, frontend.get_num_elements());
            encoder.postprocess_decode();
//            timer.stop("Decoder");

            lossless.postdecompress_data(compressed_data);

            timer.start();
            frontend.decompress(quant_inds, decData);
//            timer.stop("Prediction & Recover");
            return decData;
        }


    private:
        Frontend frontend;
        Encoder encoder;
        Lossless lossless;
    };

    template<class T, uint N, class Frontend, class Encoder, class Lossless>
    std::shared_ptr<SZAdaptiveCompressor<T, N, Frontend, Encoder, Lossless>>
    make_sz_adaptive_compressor(Frontend frontend, Encoder encoder, Lossless lossless) {
        return std::make_shared<SZAdaptiveCompressor<T, N, Frontend, Encoder, Lossless>>(frontend, encoder, lossless);
    }


}
#endif
