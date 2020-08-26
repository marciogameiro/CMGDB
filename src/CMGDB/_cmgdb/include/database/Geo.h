#ifndef GEO_H
#define GEO_H

// Base class for Geometrical object
class Geo {
public:
	virtual ~Geo ( void ) {}
  /// std stream interface
  friend std::ostream & 
  operator << ( std::ostream & out, const Geo & print_me ) {
    print_me . print( out );
    return out;
  }

  virtual std::vector<double> get_lower_bounds ( void ) const {
    std::vector<double> l_bounds = { 0.0 };
    return l_bounds;
  }

  virtual std::vector<double> get_upper_bounds ( void ) const {
    std::vector<double> u_bounds = { 0.0 };
    return u_bounds;
  }

private:
  /// derivation interface
  virtual void print ( std::ostream & ) const = 0;
};

#endif
