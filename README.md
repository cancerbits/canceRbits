To install from github:

If a `GITHUB_PAT` environment variable is already set to your access token
```
remotes::install_github(repo = 'cancerbits/canceRbits')
```

Else, pass your token
```
remotes::install_github(repo = 'cancerbits/canceRbits', auth_token = 'ghp_xxx')
```

