#ifndef __GNUPLOT_HPP
#define __GNUPLOT_HPP

#include "colsim/common.hpp"
#include "colsim/utils.hpp"

#include <vector>
#include <string>
#include <fstream>
#include <ranges>

namespace colsim
{
	template <typename T>
	requires (std::is_arithmetic_v<T>)
	class Plotter final
	{
		using range_type = float;
		using data_type = T;
		using enum_type = uint32_t;
		
	private:
		std::string _script_file_name{}; //!< name of output file
		fs::path _script_file_path{}; //!< path to script file
		std::ofstream _script_file{}; //!< filestream for handling script file
		std::string _data_file_name{}; //!< name of data file
	    fs::path _data_file_path{}; //!< path to data file
		std::ofstream _data_file{}; //!< filestream for handling data file

		static const std::string DATA_DIR;
		static const std::string SCRIPT_DIR;

		std::vector<data_type> _X{}, _y{};

		enum_type _output_file_type{}; //!< e.g. PNG, PDF
		std::string _output_file_name{};
		
		std::string _title{};
		std::string _xlabel{};
		std::string _ylabel{};

		range_type _xmin{}, _xmax{};
		range_type _ymin{}, _ymax{};

		enum_type _line_width{};
		enum_type _line_style{};
		enum_type _line_color{}; 
		enum_type _legend_loc{};
		enum_type _border_thickness{}; //!< plot outline border

		uint _plot_num{};
		std::vector<std::string> _plots{}; //!< list of created plots

	public:
		enum LineStyle : enum_type
		{
			SOLID = 0,
			DASHED
		};

		enum LegendLoc : enum_type
		{
			NONE = 0,
			TOPLEFT,
			TOPRIGHT,
			BOTTOMRIGHT,
			BOTTOMLEFT,
		};

		enum OutputFileType : enum_type
		{
			AUTO = 0, //!< guess based on extension
			PNG,
			PDF,
		};

		enum BorderThickness : enum_type
		{
			EXTRA_THIN,
			THIN,
			NORMAL,
			THICK,
			EXTRA_THICK,
		};

		enum LineColor : enum_type // must be specific here!
		{
			WHITE   = 0xFFFFFF,
			BLACK   = 0x000000,
			RED     = 0xFF0000,
			GREEN   = 0x00FF00,
			BLUE    = 0x0000FF,
			YELLOW  = 0xFFFF00,
			TEAL    = 0x00FFFF,
			MAGENTA = 0xFF00FF,
		};


	private:

		constexpr static enum_type _DEFAULT_LW = 1;
		constexpr static enum_type _DEFAULT_LC = LineColor::BLACK;
		constexpr static enum_type _DEFAULT_LS = LineStyle::SOLID;
		constexpr static enum_type _DEFAULT_LEGEND_LOC = LegendLoc::NONE;
		constexpr static enum_type _DEFAULT_BORDER_THICKNESS = BorderThickness::NORMAL;

	public:
		Plotter() = default;
		~Plotter() = default;


		
#define GEN_PLOTTER_SET_FUNC0(name, T)									\
		inline Plotter& COLSIM_JOIN(set_, name)(T const& name)			\
		{																\
			COLSIM_JOIN(_,name) = name;									\
			return *this;												\
		}
#define GEN_PLOTTER_SET_FUNC(name, T) GEN_PLOTTER_SET_FUNC0(name, T)

		GEN_PLOTTER_SET_FUNC(title, std::string)
		GEN_PLOTTER_SET_FUNC(xlabel, std::string)
		GEN_PLOTTER_SET_FUNC(ylabel, std::string)
	    GEN_PLOTTER_SET_FUNC(xmin, range_type)
		GEN_PLOTTER_SET_FUNC(xmax, range_type)
		GEN_PLOTTER_SET_FUNC(ymin, range_type)
		GEN_PLOTTER_SET_FUNC(ymax, range_type)
		GEN_PLOTTER_SET_FUNC(line_width, enum_type)
		GEN_PLOTTER_SET_FUNC(line_style, enum_type)
	    GEN_PLOTTER_SET_FUNC(legend_loc, enum_type)
		GEN_PLOTTER_SET_FUNC(border_thickness, enum_type)
		GEN_PLOTTER_SET_FUNC(line_color, enum_type)

		inline Plotter& setOutputFile(std::string const& file_name, uint output_file_type=AUTO)
		{
			if (output_file_type == AUTO) {
				std::string::size_type dot_pos = file_name.find(".");
				if (dot_pos == std::string::npos)
					log(LOG_ERROR, "Plotter::setOutputFile()", "Failed to deduce type of output file {}", file_name);

				std::string ext = file_name.substr(dot_pos + 1);

				// test if pdf
				if (ext.compare("pdf") == 0) {
					output_file_type = PDF;
				} else if (ext.compare("png") == 0) {
					output_file_type = PNG;
				} else {
					log(LOG_ERROR, "Plotter::setOutputFile()", "File extension '{}' isn't recognized.", ext);
				}
			}

			_output_file_name = file_name;
			_output_file_type = output_file_type;
		
			return *this;
		}

		Plotter& plot(std::vector<data_type> const& X, std::vector<data_type> const& y, std::string plot_title="")
		{
			if (X.size() != y.size())
				log(LOG_ERROR, "Plotter::plot()", "X size ({}) and y size ({}) differ.", X.size(), y.size());

		    _data_file_path = DATA_DIR;
			if (!fs::exists(_data_file_path)) {
				if(!fs::create_directory(_data_file_path))
					log(LOG_ERROR, "Plotter::plot()", "Failed to create output data directory '{}'.", _data_file_path.string());
			}

			std::ostringstream data_file_name{};
			data_file_name << "data" << ++_plot_num << ".dat";

			_data_file_path /= data_file_name.str();
			_data_file.open(_data_file_path.make_preferred(), std::ios_base::out);
			if (!_data_file)
				log(LOG_ERROR, "Plotter::plot()", "Failed to open file '{}'.", _data_file_path.string());

			for (std::tuple<data_type,data_type> x : std::ranges::zip_view{X, y})
				_data_file << std::get<0>(x) << '\t' << std::get<1>(x) << '\n';
		    
			_X = X;
			_y = y;
			_title = plot_title;

			_data_file.close();

			return *this;
		}

		inline Plotter& save()
		{
			_script_file_path = SCRIPT_DIR;
			if (!fs::exists(_script_file_path)) {
				if (!fs::create_directory(_script_file_path))
					log(LOG_ERROR, "Plotter::plot()", "Failed to create script output directory.");
			}

			_script_file_name = "script.gplt";
			_script_file_path /= _script_file_name;
			_script_file.open(_script_file_path);
			if (!_script_file)
				log(LOG_ERROR, "Plotter::save()", "Failed to open the script file '{}'.", _script_file_path.string());

			std::string terminal_type{};
			switch (_output_file_type) {
				case PNG: terminal_type = "pngcairo"; break;
				case PDF: terminal_type = "pdfcairo"; break;
				case AUTO:
				default: log(LOG_ERROR, "Plotter::save()", "Invalid output file/terminal type.");
			}

			_script_file << "set terminal " << terminal_type << " enhanced notransparent\n";
			_script_file << "set output '" << _output_file_name << "'\n";
			
			// I could throw an error if the user doesn't specify a valid range,
			// but instead I'll just iterate through the provided lists
			// and pick out the actual min and max +- 3%
			if ((_xmin == 0 and _xmax == 0) || _xmin == _xmax) {
				range_type xmin = *std::min_element(_X.begin(), _X.end());
				range_type xmax = *std::max_element(_X.begin(), _X.end());

				range_type three_percent = (xmax-xmin)*static_cast<range_type>(0.03);
				// three_percent = 0.0;
				
				_xmin = xmin - three_percent;
				_xmax = xmax + three_percent;
			}
			if ((_ymin == 0 and _ymax == 0) || _ymin == _ymax) {
				range_type ymin = *std::min_element(_y.begin(), _y.end());
				range_type ymax = *std::max_element(_y.begin(), _y.end());

				range_type three_percent = (ymax-ymin)*static_cast<range_type>(0.03);
				// three_percent = 0.0;
				
				_ymin = ymin - three_percent;
				_ymax = ymax + three_percent;
			}
			_script_file << _genLine_xrange() << '\n';
			_script_file << _genLine_yrange() << '\n';
			_script_file << "set xlabel '" << _xlabel << "'\n";
			_script_file << "set ylabel '" << _ylabel << "'\n";
			_script_file << "set title '" << _title << "'\n";

			_script_file << "plot '" << _data_file_path.string() << "' with lines title '" << _title << "'\n";
			_script_file.close();

			std::ostringstream command{};
			command << "gnuplot " << _script_file_path.string();
			log(LOG_INFO, "Plotter::save()", "the command is {}", command.str());

			if (system(command.str().c_str()) != 0)
				log(LOG_ERROR, "Plotter::save()", "Failed to create child shell to call gnuplot, or gnuplot call failed.");

			return *this;
		}


	private:
		inline std::string _genLine_yrange()
		{
			if (_ymin >= _ymax)
				log(LOG_ERROR, "Plotter::_genLine_yrange()", "minimum y value ({}) is equal to or larger than maximum y value ({})", _ymin, _ymax);
		
			std::ostringstream ss{};
			ss << "set yrange [" << _ymin << ":" << _ymax << "]";
			return ss.str();
		}
		
		inline std::string _genLine_xrange()
		{
			if (_xmin >= _xmax)
				log(LOG_ERROR, "Plotter::_genLine_xrange()", "minimum x value ({}) is equal to or larger than maximum x value ({})", _xmin, _xmax);
		
			std::ostringstream ss{};
			ss << "set yrange [" << _xmin << ":" << _xmax << "]";
			return ss.str();
		}
	};


} // namespace colsim






#endif // endif __GNUPLOT_HPP
