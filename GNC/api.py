import requests
import sys
def getCurrentWeather(lat,lon):
    base_url = "https://api.openweathermap.org/data/2.5/weather"
    key = "a71b6c2361abf47745ec9cd4bb69d470" 
    units = "imperial"
    key_param = f"appid={key}"
    lat_param = f"lat={lat}"
    lon_param = f"lon={lon}"
    units_param= f"units={units}"
    query_string = f"{units_param}&{lat_param}&{lon_param}&{key_param}"
    full_url = f"{base_url}?{query_string}"
    response = requests.get(full_url)
    data = response.json()
    country = data['sys']['country']
    city = data['name']
    air_pressure_sea_lvl = data['main']['sea_level']
    wind_speed = data['wind']['speed']
    wind_direction = data['wind']['deg']
    weather = '['
    for i in range(len(data['weather'])):
        if i != len(data['weather'])-1:
            weather += data['weather'][i]['main']+"-"+data['weather'][i]['description']+"|"
        else:
            weather += data['weather'][i]['main']+"-"+data['weather'][i]['description']
    weather +=']'

    print(f'{country},{city},{weather},{air_pressure_sea_lvl},{wind_speed},{wind_direction}')

lat = float(sys.argv[1])
lon = float(sys.argv[2])

getCurrentWeather(lat,lon)