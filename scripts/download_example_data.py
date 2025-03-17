import requests

response = requests.get('https://byu.box.com/shared/static/seu2e75iw7ihnxqnsp8gi29xmo9pqynr.zip', stream=True)
with open('mass_spec.zip', 'wb') as f:
    for chunk in response.iter_content(chunk_size=8192):
        f.write(chunk)
