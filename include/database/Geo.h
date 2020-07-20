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
private:
  /// derivation interface
  virtual void print ( std::ostream & ) const = 0;
};

#endif
